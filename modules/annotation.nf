/*
========================================================================================
    ANNOTATION MODULE
========================================================================================
    This module handles the annotation of filtered variants:
    - Variant Effect Predictor (VEP) for functional annotation
    - Integration with clinical databases (ClinVar, COSMIC, etc.)
    - Custom annotation with gene panels and disease databases
    - Structural variant annotation
*/

// Module parameters with defaults
params.vep_cache = null
params.vep_plugins = "Mastermind,CADD,SpliceAI,dbNSFP"
params.clinvar = null
params.cosmic = null
params.gnomad = null
params.custom_annotations = null
params.cache_version = "108"
params.species = "homo_sapiens"
params.assembly = "GRCh38"
params.vep_fields = "Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,CANONICAL,ClinVar,Mastermind_URL,CADD_PHRED,CADD_RAW,SpliceAI_pred_DP_AG,SpliceAI_pred_DP_AL,SpliceAI_pred_DP_DG,SpliceAI_pred_DP_DL,SpliceAI_pred_DS_AG,SpliceAI_pred_DS_AL,SpliceAI_pred_DS_DG,SpliceAI_pred_DS_DL,gnomAD_AF,gnomAD_AF_popmax"

/*
 * PROCESS: Run VEP for functional annotation
 */
process VEP {
    tag "$meta.id"
    label 'process_medium'
    container 'ensemblorg/ensembl-vep:release_108.1'
    
    publishDir "${params.output}/variants/annotated", mode: 'copy',
        saveAs: { filename -> 
            if (filename.endsWith('.summary.html')) "vep_summaries/$filename"
            else if (filename.endsWith('.json')) "json/$filename"
            else "$filename"
        }
    
    input:
    tuple val(meta), path(vcf), path(tbi)
    path vep_cache
    path clinvar
    path cosmic
    path gnomad
    path custom_annotations
    
    output:
    tuple val(meta), path("*.ann.vcf.gz"), emit: vcf
    tuple val(meta), path("*.ann.vcf.gz.tbi"), emit: tbi
    tuple val(meta), path("*.summary.html"), emit: html
    tuple val(meta), path("*.ann.json"), emit: json
    path "versions.yml", emit: versions
    
    script:
    def prefix = meta.id
    def cache_option = vep_cache ? "--cache --dir_cache ${vep_cache}" : "--database"
    def plugin_option = params.vep_plugins ? "--plugin ${params.vep_plugins}" : ""
    def clinvar_option = clinvar ? "--custom ${clinvar},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN" : ""
    def cosmic_option = cosmic ? "--custom ${cosmic},COSMIC,vcf,exact,0,COSV,GENE,CDS,AA" : ""
    def gnomad_option = gnomad ? "--custom ${gnomad},gnomAD,vcf,exact,0,AF,AF_popmax" : ""
    def custom_option = custom_annotations ? "--custom ${custom_annotations.join(',')}" : ""
    
    """
    vep \\
        --input_file ${vcf} \\
        --output_file ${prefix}.ann.vcf.gz \\
        --vcf \\
        --compress_output bgzip \\
        --species ${params.species} \\
        --assembly ${params.assembly} \\
        --cache_version ${params.cache_version} \\
        --fields "${params.vep_fields}" \\
        --json \\
        --force_overwrite \\
        --stats_file ${prefix}.summary.html \\
        ${cache_option} \\
        ${plugin_option} \\
        ${clinvar_option} \\
        ${cosmic_option} \\
        ${gnomad_option} \\
        ${custom_option} \\
        --fork ${task.cpus}
    
    # Create json output
    bcftools view ${prefix}.ann.vcf.gz | grep -v "^#" | jq -R -s 'split("\\n") | map(select(length > 0)) | map(split("\\t")) | map({CHROM: .[0], POS: .[1], ID: .[2], REF: .[3], ALT: .[4], QUAL: .[5], FILTER: .[6], INFO: .[7], FORMAT: .[8]})' > ${prefix}.ann.json
    
    # Index the VCF
    tabix -p vcf ${prefix}.ann.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$(vep --help 2>&1 | grep "ensembl-vep" | sed 's/.*ensembl-vep : //;s/ .*\$//')
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools //g')
        tabix: \$(tabix --version | head -n1 | sed 's/tabix (htslib) //g')
    END_VERSIONS
    """
}

/*
 * PROCESS: Post-process VEP output to add clinical significance
 */
process CLINICAL_ANNOTATION {
    tag "$meta.id"
    label 'process_low'
    container 'python:3.9'
    
    publishDir "${params.output}/variants/annotated", mode: 'copy'
    
    input:
    tuple val(meta), path(vcf), path(tbi)
    path gene_panel
    
    output:
    tuple val(meta), path("*.clinical.vcf.gz"), emit: vcf
    tuple val(meta), path("*.clinical.vcf.gz.tbi"), emit: tbi
    tuple val(meta), path("*.clinical.json"), emit: json
    path "versions.yml", emit: versions
    
    script:
    def prefix = meta.id
    def gene_panel_arg = gene_panel ? "--gene_panel $gene_panel" : ""
    
    """
    #!/usr/bin/env python3
    
    import gzip
    import json
    import sys
    import os
    import re
    from pathlib import Path
    
    # Sample Python implementation that would:
    # 1. Parse input VCF
    # 2. Add clinical annotations based on gene panel data
    # 3. Output annotated VCF and JSON
    
    # In a real implementation, this would be a more sophisticated module
    # that integrates with various clinical knowledge bases
    
    # For demonstration purposes, just copy the input and add minimal annotation
    with gzip.open("${vcf}", 'rt') as in_vcf, gzip.open("${prefix}.clinical.vcf.gz", 'wt') as out_vcf:
        for line in in_vcf:
            if line.startswith('#'):
                if line.startswith('##INFO'):
                    out_vcf.write(line)
                    if '##INFO=<ID=CSQ' in line:
                        # Add custom clinical annotation fields
                        out_vcf.write('##INFO=<ID=CLINICAL_SIG,Number=.,Type=String,Description="Clinical significance">')
                        out_vcf.write('##INFO=<ID=IN_PANEL,Number=0,Type=Flag,Description="Variant in gene panel">')
                else:
                    out_vcf.write(line)
            else:
                # Process variant line (add clinical annotations)
                parts = line.strip().split('\\t')
                chrom, pos, id, ref, alt, qual, filter, info = parts[:8]
                
                # Add clinical annotations (this would be more sophisticated in real implementation)
                info += ";CLINICAL_SIG=Uncertain_significance;IN_PANEL"
                
                parts[7] = info
                out_vcf.write('\\t'.join(parts) + '\\n')
    
    # Create JSON output with clinical annotation
    with gzip.open("${vcf}", 'rt') as in_vcf:
        variants = []
        for line in in_vcf:
            if not line.startswith('#'):
                parts = line.strip().split('\\t')
                if len(parts) >= 8:
                    chrom, pos, id, ref, alt = parts[:5]
                    info = parts[7]
                    
                    # Extract gene from CSQ field
                    gene = "Unknown"
                    csq_match = re.search(r'CSQ=([^;]+)', info)
                    if csq_match:
                        csq_parts = csq_match.group(1).split('|')
                        if len(csq_parts) > 3:
                            gene = csq_parts[3]
                    
                    # Create variant record
                    variant = {
                        "chrom": chrom,
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                        "gene": gene,
                        "clinical_significance": "Uncertain_significance",
                        "in_panel": True
                    }
                    variants.append(variant)
        
        # Write to JSON
        with open("${prefix}.clinical.json", 'w') as out_json:
            json.dump({"sample_id": "${meta.id}", "variants": variants}, out_json, indent=2)
    
    # Index output VCF
    os.system("tabix -p vcf ${prefix}.clinical.vcf.gz")
    
    with open("versions.yml", "w") as versions_file:
        versions_file.write(\"""${task.process}":
        python: \$(python --version | sed 's/Python //g')
        tabix: \$(tabix --version | head -n1 | sed 's/tabix (htslib) //g')
    \""")
    """
}

/*
 * PROCESS: Generate annotated variant summary
 */
process VARIANT_SUMMARY {
    tag "$meta.id"
    label 'process_low'
    container 'python:3.9'
    
    publishDir "${params.output}/variants/annotated", mode: 'copy'
    
    input:
    tuple val(meta), path(vcf), path(tbi)
    
    output:
    tuple val(meta), path("*.variants_summary.csv"), emit: csv
    tuple val(meta), path("*.variants_summary.json"), emit: json
    path "versions.yml", emit: versions
    
    script:
    def prefix = meta.id
    
    """
    #!/usr/bin/env python3
    
    import gzip
    import json
    import csv
    from pathlib import Path
    
    # In a real implementation, this would:
    # 1. Parse the annotated VCF
    # 2. Extract key variant information
    # 3. Generate summary in CSV and JSON formats
    
    variants = []
    
    # Parse VCF to extract variant information
    with gzip.open("${vcf}", 'rt') as vcf_file:
        for line in vcf_file:
            if not line.startswith('#'):
                parts = line.strip().split('\\t')
                if len(parts) >= 8:
                    chrom = parts[0]
                    pos = parts[1]
                    id = parts[2]
                    ref = parts[3]
                    alt = parts[4]
                    qual = parts[5]
                    filter = parts[6]
                    info = parts[7]
                    
                    # Extract annotations from info field
                    gene = "Unknown"
                    consequence = "Unknown"
                    clinical_sig = "Unknown"
                    
                    if "CSQ=" in info:
                        csq_parts = info.split("CSQ=")[1].split(";")[0].split("|")
                        if len(csq_parts) > 3:
                            consequence = csq_parts[0] 
                            gene = csq_parts[3]
                    
                    if "CLINICAL_SIG=" in info:
                        clinical_sig = info.split("CLINICAL_SIG=")[1].split(";")[0]
                    
                    # Create variant record
                    variant = {
                        "chrom": chrom,
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                        "gene": gene,
                        "consequence": consequence,
                        "clinical_significance": clinical_sig,
                        "filter": filter
                    }
                    variants.append(variant)
    
    # Write variants to CSV
    with open("${prefix}.variants_summary.csv", 'w', newline='') as csv_file:
        fieldnames = ["chrom", "pos", "ref", "alt", "gene", "consequence", "clinical_significance", "filter"]
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(variants)
    
    # Write variants to JSON
    with open("${prefix}.variants_summary.json", 'w') as json_file:
        json.dump({"sample_id": "${meta.id}", "variants": variants}, json_file, indent=2)
    
    with open("versions.yml", "w") as versions_file:
        versions_file.write(\"""${task.process}":
        python: \$(python --version | sed 's/Python //g')
    \""")
    """
}

// Main workflow for annotation
workflow ANNOTATION {
    take:
    vcf                // [meta, vcf, tbi]
    vep_cache          // VEP cache directory
    clinvar            // ClinVar VCF file
    cosmic             // COSMIC VCF file
    gnomad             // gnomAD VCF file
    custom_annotations // List of custom annotation files
    gene_panel         // Gene panel JSON file
    
    main:
    ch_versions = Channel.empty()
    
    // Run VEP
    VEP(
        vcf,
        vep_cache,
        clinvar,
        cosmic,
        gnomad,
        custom_annotations
    )
    ch_versions = ch_versions.mix(VEP.out.versions)
    
    // Add clinical annotations
    ch_vep_vcf_tbi = VEP.out.vcf.join(VEP.out.tbi, by: [0])
    CLINICAL_ANNOTATION(
        ch_vep_vcf_tbi,
        gene_panel
    )
    ch_versions = ch_versions.mix(CLINICAL_ANNOTATION.out.versions)
    
    // Generate variant summary
    ch_clinical_vcf_tbi = CLINICAL_ANNOTATION.out.vcf.join(CLINICAL_ANNOTATION.out.tbi, by: [0])
    VARIANT_SUMMARY(
        ch_clinical_vcf_tbi
    )
    ch_versions = ch_versions.mix(VARIANT_SUMMARY.out.versions)
    
    emit:
    vcf = CLINICAL_ANNOTATION.out.vcf
    tbi = CLINICAL_ANNOTATION.out.tbi
    json = CLINICAL_ANNOTATION.out.json
    summary_csv = VARIANT_SUMMARY.out.csv
    summary_json = VARIANT_SUMMARY.out.json
    vep_html = VEP.out.html
    versions = ch_versions
} 