#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Report Generator for Clinical Pipeline

This module generates clinical reports from variant data using templates.
"""

import os
import sys
import json
import logging
import argparse
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Union

import jinja2
import jsonschema
import markdown
import pandas as pd
from weasyprint import HTML, CSS

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger("REPORT_GENERATOR")

class ClinicalReportGenerator:
    """Generator for clinical genomic reports."""
    
    def __init__(
        self, 
        template_dir: Union[str, Path], 
        schema_path: Optional[Union[str, Path]] = None,
        anonymize: bool = False
    ):
        """
        Initialize the report generator.
        
        Args:
            template_dir: Directory containing report templates
            schema_path: Path to JSON schema for report data validation
            anonymize: Whether to generate anonymized reports
        """
        self.template_dir = Path(template_dir)
        self.anonymize = anonymize
        
        # Set up Jinja2 environment
        self.jinja_env = jinja2.Environment(
            loader=jinja2.FileSystemLoader(str(self.template_dir)),
            autoescape=jinja2.select_autoescape(['html', 'xml'])
        )
        
        # Load schema if provided
        self.schema = None
        if schema_path:
            with open(schema_path, 'r') as f:
                self.schema = json.load(f)
    
    def load_data(self, data_path: Union[str, Path]) -> Dict:
        """
        Load and validate report data from JSON file.
        
        Args:
            data_path: Path to JSON data file
            
        Returns:
            Dictionary containing report data
            
        Raises:
            jsonschema.exceptions.ValidationError: If data doesn't match schema
        """
        with open(data_path, 'r') as f:
            data = json.load(f)
        
        # Validate data against schema if available
        if self.schema:
            jsonschema.validate(instance=data, schema=self.schema)
        
        return data
    
    def load_variants_from_tsv(self, tsv_path: Union[str, Path]) -> List[Dict]:
        """
        Load variant data from TSV file.
        
        Args:
            tsv_path: Path to TSV file containing variant data
            
        Returns:
            List of dictionaries containing variant data
        """
        df = pd.read_csv(tsv_path, sep='\t')
        
        variants = []
        for _, row in df.iterrows():
            # Map database fields to report fields
            variant = {
                'gene': row.get('GENE', '.'),
                'hgvs': f"{row.get('HGVSc', '.')} ({row.get('HGVSp', '.')})",
                'zygosity': self._determine_zygosity(row),
                'transcript': row.get('Feature', '.'),
                'significance': self._get_clinical_significance(row),
                'classification': self._map_classification(row.get('CLINVAR', '.')),
                'significance_class': self._get_significance_class(row.get('CLINVAR', '.'))
            }
            variants.append(variant)
        
        return variants
    
    def _determine_zygosity(self, variant_row: pd.Series) -> str:
        """Determine zygosity from variant data."""
        af = variant_row.get('AF', 0)
        if isinstance(af, str) and af == '.':
            return 'Unknown'
        
        try:
            af = float(af)
            if af <= 0.6:
                return 'Heterozygous'
            elif af > 0.8:
                return 'Homozygous'
            else:
                return 'Mosaic'
        except (ValueError, TypeError):
            return 'Unknown'
    
    def _get_clinical_significance(self, variant_row: pd.Series) -> str:
        """Generate clinical significance description from variant data."""
        impact = variant_row.get('IMPACT', '')
        consequence = variant_row.get('CONSEQUENCE', '')
        clinvar = variant_row.get('CLINVAR', '')
        
        if clinvar and clinvar != '.':
            return f"This variant is classified as {clinvar.lower()} according to ClinVar."
        elif impact == 'HIGH':
            return f"This variant has a high impact ({consequence}) on the protein function."
        elif impact == 'MODERATE':
            return f"This variant has a moderate impact ({consequence}) on the protein function."
        else:
            return "The clinical significance of this variant is uncertain."
    
    def _map_classification(self, clinvar: str) -> str:
        """Map ClinVar classification to report classification."""
        mapping = {
            'pathogenic': 'Pathogenic',
            'likely_pathogenic': 'Likely Pathogenic',
            'uncertain_significance': 'Uncertain Significance',
            'likely_benign': 'Likely Benign',
            'benign': 'Benign'
        }
        
        # Try direct match
        if clinvar.lower() in mapping:
            return mapping[clinvar.lower()]
        
        # Try substring match
        for key, value in mapping.items():
            if key.replace('_', ' ') in clinvar.lower():
                return value
        
        return 'Uncertain Significance'
    
    def _get_significance_class(self, classification: str) -> str:
        """Convert classification to CSS class name."""
        classification = classification.lower()
        
        if 'pathogenic' in classification and 'likely' not in classification:
            return 'pathogenic'
        elif 'likely pathogenic' in classification or 'likely_pathogenic' in classification:
            return 'likely-pathogenic'
        elif 'benign' in classification and 'likely' not in classification:
            return 'benign'
        elif 'likely benign' in classification or 'likely_benign' in classification:
            return 'likely-benign'
        else:
            return 'vus'
    
    def generate_html_report(
        self, 
        report_data: Dict, 
        output_path: Union[str, Path]
    ) -> Path:
        """
        Generate HTML report from template and data.
        
        Args:
            report_data: Dictionary containing report data
            output_path: Path where HTML report will be saved
            
        Returns:
            Path to the generated HTML report
        """
        template_name = 'anonymous_report_template.html' if self.anonymize else 'report_template.html'
        template = self.jinja_env.get_template(template_name)
        
        # Render the template with report data
        html_content = template.render(**report_data)
        
        # Write HTML file
        output_path = Path(output_path)
        with open(output_path, 'w') as f:
            f.write(html_content)
        
        logger.info(f"HTML report generated at {output_path}")
        return output_path
    
    def generate_pdf_report(
        self, 
        html_path: Union[str, Path],
        output_path: Union[str, Path],
        css_path: Optional[Union[str, Path]] = None
    ) -> Path:
        """
        Generate PDF report from HTML report.
        
        Args:
            html_path: Path to HTML report
            output_path: Path where PDF report will be saved
            css_path: Path to CSS file for styling (optional)
            
        Returns:
            Path to the generated PDF report
        """
        html_path = Path(html_path)
        output_path = Path(output_path)
        
        # Load HTML file
        html = HTML(filename=str(html_path))
        
        # Apply CSS if provided
        stylesheets = []
        if css_path:
            stylesheets.append(CSS(filename=str(css_path)))
        
        # Generate PDF
        html.write_pdf(output_path, stylesheets=stylesheets)
        
        logger.info(f"PDF report generated at {output_path}")
        return output_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Clinical Report Generator")
    parser.add_argument('--template_dir', required=True, help='Directory containing report templates')
    parser.add_argument('--data', required=True, help='Path to report data JSON or variants TSV')
    parser.add_argument('--output_html', required=True, help='Output HTML report path')
    parser.add_argument('--output_pdf', help='Output PDF report path (optional)')
    parser.add_argument('--schema', help='JSON schema for data validation (optional)')
    parser.add_argument('--css', help='CSS file for PDF styling (optional)')
    parser.add_argument('--anonymize', action='store_true', help='Generate anonymized report')
    
    args = parser.parse_args()
    
    try:
        # Initialize report generator
        report_generator = ClinicalReportGenerator(
            template_dir=args.template_dir,
            schema_path=args.schema,
            anonymize=args.anonymize
        )
        
        # Determine if input is JSON or TSV
        input_path = Path(args.data)
        if input_path.suffix.lower() == '.json':
            # Load data from JSON
            report_data = report_generator.load_data(input_path)
        elif input_path.suffix.lower() in ('.tsv', '.txt'):
            # Load data from TSV
            # For simplicity, we'll create a minimal report data structure
            variants = report_generator.load_variants_from_tsv(input_path)
            report_data = {
                'report_id': f"GR-{datetime.now().strftime('%Y%m%d%H%M%S')}",
                'report_title': 'Genomic Variant Report',
                'report_date': datetime.now().strftime('%Y-%m-%d'),
                'patient_id': 'ANONYMOUS' if args.anonymize else 'PATIENT',
                'sample_id': Path(input_path).stem,
                'results_summary': f"Found {len(variants)} variants of interest.",
                'variants': variants,
                'clinical_interpretation': "Please interpret these findings in context of clinical presentation.",
                'result_class': 'positive' if any(v['significance_class'] in ['pathogenic', 'likely-pathogenic'] for v in variants) else 'inconclusive'
            }
        else:
            raise ValueError(f"Unsupported input file format: {input_path.suffix}")
        
        # Generate HTML report
        html_path = report_generator.generate_html_report(
            report_data=report_data,
            output_path=args.output_html
        )
        
        # Generate PDF report if requested
        if args.output_pdf:
            report_generator.generate_pdf_report(
                html_path=html_path,
                output_path=args.output_pdf,
                css_path=args.css
            )
        
        logger.info("Report generation completed successfully")
        
    except Exception as e:
        logger.error(f"Error generating report: {str(e)}")
        sys.exit(1) 