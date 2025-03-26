#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Variant Filtering Script for Clinical Pipeline

This script implements clinical-grade variant filtering for the bioinformatics pipeline,
applying frequency, quality, and impact filters to identify clinically significant variants.
"""

import os
import sys
import logging
import argparse
from pathlib import Path
from typing import Dict, List, Optional, Set, Union

import pysam
import pandas as pd
import numpy as np
from tqdm import tqdm

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger("VARIANT_FILTER")

def read_panel_genes(panel_file: Path) -> Set[str]:
    """
    Read gene panel from file (one gene per line).

    Args:
        panel_file: Path to gene panel file

    Returns:
        Set of gene names
    """
    panel_genes = set()
    with open(panel_file, 'r') as f:
        for line in f:
            gene = line.strip()
            if gene and not gene.startswith('#'):
                panel_genes.add(gene)
    
    logger.info(f"Loaded {len(panel_genes)} genes from panel file")
    return panel_genes

def parse_vep_annotation(info_field: Dict) -> Dict:
    """
    Parse VEP annotation from VCF INFO field.

    Args:
        info_field: INFO field from VCF

    Returns:
        Dictionary with parsed VEP annotations
    """
    if 'CSQ' not in info_field:
        return {}
    
    # Get VEP annotations and field names
    csq_data = info_field['CSQ']
    if not csq_data:
        return {}
    
    # VEP format is in the INFO header
    # Usually: Consequence|IMPACT|SYMBOL|Gene|...
    csq_fields = ['Consequence', 'IMPACT', 'SYMBOL', 'Gene', 'Feature_type', 
                 'Feature', 'BIOTYPE', 'EXON', 'INTRON', 'HGVSc', 'HGVSp', 
                 'cDNA_position', 'CDS_position', 'Protein_position', 
                 'Amino_acids', 'Codons', 'STRAND', 'SIFT', 'PolyPhen',
                 'gnomAD_AF', 'CLIN_SIG', 'PUBMED']
    
    # Parse the most severe consequence (first in list)
    annotations = csq_data.split(',')[0].split('|')
    
    # Create dictionary of annotations
    vep_dict = {csq_fields[i]: annotations[i] for i in range(min(len(csq_fields), len(annotations)))}
    
    return vep_dict

def filter_variants(
    vcf_path: Path,
    output_vcf: Path,
    output_tsv: Path,
    panel_genes: Optional[Set[str]] = None,
    maf_cutoff: float = 0.01,
    min_depth: int = 20,
    min_quality: float = 20.0
):
    """
    Filter VCF file based on clinical criteria.

    Args:
        vcf_path: Path to input VCF file
        output_vcf: Path to output filtered VCF file
        output_tsv: Path to output TSV file with filtered variants
        panel_genes: Set of genes to include (optional)
        maf_cutoff: Minor allele frequency cutoff
        min_depth: Minimum read depth
        min_quality: Minimum variant quality
    """
    # Open input VCF
    vcf_in = pysam.VariantFile(str(vcf_path))
    
    # Create output VCF with same header
    vcf_out = pysam.VariantFile(str(output_vcf), 'w', header=vcf_in.header)
    
    # Prepare TSV output
    tsv_columns = ['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER', 
                   'GENE', 'IMPACT', 'CONSEQUENCE', 'HGVSc', 'HGVSp',
                   'gnomAD_AF', 'DP', 'AF', 'CLINVAR']
    tsv_data = []
    
    # Process variants
    total_variants = 0
    filtered_variants = 0
    
    logger.info(f"Filtering variants in {vcf_path}")
    for record in tqdm(vcf_in.fetch(), desc="Processing variants"):
        total_variants += 1
        
        # Basic quality filters
        if record.qual is not None and record.qual < min_quality:
            continue
        
        # Extract depth and allele frequency for each sample
        for sample in record.samples:
            if 'DP' in record.samples[sample] and record.samples[sample]['DP'] < min_depth:
                continue
        
        # Extract VEP annotations
        vep_annotation = parse_vep_annotation(record.info)
        
        # Panel gene filter
        if panel_genes and vep_annotation.get('SYMBOL') not in panel_genes:
            continue
        
        # Frequency filter
        gnomad_af = vep_annotation.get('gnomAD_AF', '.')
        if gnomad_af != '.' and float(gnomad_af) > maf_cutoff:
            continue
        
        # Impact filter - keep high/moderate impact variants
        impact = vep_annotation.get('IMPACT', '')
        if impact not in ['HIGH', 'MODERATE']:
            continue
        
        # Passed all filters
        filtered_variants += 1
        
        # Add to output VCF
        vcf_out.write(record)
        
        # Add to TSV data
        for alt in record.alts:
            tsv_row = {
                'CHROM': record.chrom,
                'POS': record.pos,
                'REF': record.ref,
                'ALT': alt,
                'QUAL': record.qual,
                'FILTER': ','.join(record.filter.keys()) if record.filter.keys() else 'PASS',
                'GENE': vep_annotation.get('SYMBOL', '.'),
                'IMPACT': vep_annotation.get('IMPACT', '.'),
                'CONSEQUENCE': vep_annotation.get('Consequence', '.'),
                'HGVSc': vep_annotation.get('HGVSc', '.'),
                'HGVSp': vep_annotation.get('HGVSp', '.'),
                'gnomAD_AF': vep_annotation.get('gnomAD_AF', '.'),
                'DP': record.info.get('DP', '.'),
                'AF': record.info.get('AF', ['.'])[0],
                'CLINVAR': vep_annotation.get('CLIN_SIG', '.')
            }
            tsv_data.append(tsv_row)
    
    # Close VCF files
    vcf_in.close()
    vcf_out.close()
    
    # Create index for output VCF
    pysam.tabix_index(str(output_vcf), preset='vcf', force=True)
    
    # Write TSV output
    if tsv_data:
        df = pd.DataFrame(tsv_data)
        df.to_csv(output_tsv, sep='\t', index=False)
    else:
        # Create empty TSV with headers
        with open(output_tsv, 'w') as f:
            f.write('\t'.join(tsv_columns) + '\n')
    
    logger.info(f"Filtering complete: {filtered_variants} / {total_variants} variants passed filters")
    logger.info(f"Results written to {output_vcf} and {output_tsv}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Clinical variant filtering")
    parser.add_argument('--vcf', required=True, help='Input VCF file')
    parser.add_argument('--output', required=True, help='Output VCF file')
    parser.add_argument('--tsv', required=True, help='Output TSV file')
    parser.add_argument('--panel', help='Gene panel file (one gene per line)')
    parser.add_argument('--maf_cutoff', type=float, default=0.01, help='MAF cutoff (default: 0.01)')
    parser.add_argument('--min_depth', type=int, default=20, help='Minimum read depth (default: 20)')
    parser.add_argument('--min_quality', type=float, default=20.0, help='Minimum variant quality (default: 20.0)')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads to use')
    
    args = parser.parse_args()
    
    try:
        # Read gene panel if provided
        panel_genes = None
        if args.panel:
            panel_genes = read_panel_genes(Path(args.panel))
        
        # Filter variants
        filter_variants(
            vcf_path=Path(args.vcf),
            output_vcf=Path(args.output),
            output_tsv=Path(args.tsv),
            panel_genes=panel_genes,
            maf_cutoff=args.maf_cutoff,
            min_depth=args.min_depth,
            min_quality=args.min_quality
        )
        
    except Exception as e:
        logger.error(f"Error in variant filtering: {str(e)}")
        sys.exit(1) 