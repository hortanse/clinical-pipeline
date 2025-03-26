#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Variant Classifier for Clinical Pipeline

This module implements ACMG/AMP guideline-based variant classification logic.
"""

import os
import sys
import json
import logging
import argparse
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union

import pandas as pd

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger("VARIANT_CLASSIFIER")

class ACMGClassifier:
    """
    ACMG/AMP guideline-based variant classifier.
    
    Implements the 2015 ACMG/AMP guidelines for variant classification:
    Richards et al. Genetics in Medicine (2015) 17, 405–424
    """
    
    # ACMG criteria weights
    CRITERIA_WEIGHTS = {
        # Pathogenic - Very Strong
        'PVS1': 'Very Strong',
        
        # Pathogenic - Strong
        'PS1': 'Strong',
        'PS2': 'Strong',
        'PS3': 'Strong',
        'PS4': 'Strong',
        
        # Pathogenic - Moderate
        'PM1': 'Moderate',
        'PM2': 'Moderate',
        'PM3': 'Moderate',
        'PM4': 'Moderate',
        'PM5': 'Moderate',
        'PM6': 'Moderate',
        
        # Pathogenic - Supporting
        'PP1': 'Supporting',
        'PP2': 'Supporting',
        'PP3': 'Supporting',
        'PP4': 'Supporting',
        'PP5': 'Supporting',
        
        # Benign - Stand-Alone
        'BA1': 'Stand-Alone',
        
        # Benign - Strong
        'BS1': 'Strong',
        'BS2': 'Strong',
        'BS3': 'Strong',
        'BS4': 'Strong',
        
        # Benign - Supporting
        'BP1': 'Supporting',
        'BP2': 'Supporting',
        'BP3': 'Supporting',
        'BP4': 'Supporting',
        'BP5': 'Supporting',
        'BP6': 'Supporting',
        'BP7': 'Supporting'
    }
    
    def __init__(self, rules_path: Optional[Union[str, Path]] = None):
        """
        Initialize the ACMG classifier.
        
        Args:
            rules_path: Path to custom ACMG rules JSON file (optional)
        """
        # Default classification rules based on ACMG/AMP guidelines
        self.classification_rules = {
            'pathogenic': [
                # Criteria for Pathogenic classification
                ['PVS1', 'PS1|PS2|PS3|PS4'],
                ['PVS1', 'PM1|PM2|PM3|PM4|PM5|PM6', 'PM1|PM2|PM3|PM4|PM5|PM6'],
                ['PVS1', 'PM1|PM2|PM3|PM4|PM5|PM6', 'PP1|PP2|PP3|PP4|PP5'],
                ['PS1|PS2|PS3|PS4', 'PS1|PS2|PS3|PS4'],
                ['PS1|PS2|PS3|PS4', 'PM1|PM2|PM3|PM4|PM5|PM6', 'PM1|PM2|PM3|PM4|PM5|PM6'],
                ['PS1|PS2|PS3|PS4', 'PM1|PM2|PM3|PM4|PM5|PM6', 'PP1|PP2|PP3|PP4|PP5', 'PP1|PP2|PP3|PP4|PP5']
            ],
            'likely_pathogenic': [
                # Criteria for Likely Pathogenic classification
                ['PVS1', 'PM1|PM2|PM3|PM4|PM5|PM6'],
                ['PVS1', 'PP1|PP2|PP3|PP4|PP5', 'PP1|PP2|PP3|PP4|PP5'],
                ['PS1|PS2|PS3|PS4', 'PM1|PM2|PM3|PM4|PM5|PM6'],
                ['PS1|PS2|PS3|PS4', 'PP1|PP2|PP3|PP4|PP5', 'PP1|PP2|PP3|PP4|PP5'],
                ['PM1|PM2|PM3|PM4|PM5|PM6', 'PM1|PM2|PM3|PM4|PM5|PM6', 'PM1|PM2|PM3|PM4|PM5|PM6'],
                ['PM1|PM2|PM3|PM4|PM5|PM6', 'PM1|PM2|PM3|PM4|PM5|PM6', 'PP1|PP2|PP3|PP4|PP5', 'PP1|PP2|PP3|PP4|PP5']
            ],
            'benign': [
                # Criteria for Benign classification
                ['BA1'],
                ['BS1|BS2|BS3|BS4', 'BS1|BS2|BS3|BS4']
            ],
            'likely_benign': [
                # Criteria for Likely Benign classification
                ['BS1|BS2|BS3|BS4', 'BP1|BP2|BP3|BP4|BP5|BP6|BP7'],
                ['BP1|BP2|BP3|BP4|BP5|BP6|BP7', 'BP1|BP2|BP3|BP4|BP5|BP6|BP7']
            ]
        }
        
        # Load custom rules if provided
        if rules_path:
            with open(rules_path, 'r') as f:
                custom_rules = json.load(f)
                self.classification_rules.update(custom_rules)
    
    def classify_variant(self, acmg_criteria: Dict[str, bool]) -> Tuple[str, List[str]]:
        """
        Classify a variant based on ACMG criteria.
        
        Args:
            acmg_criteria: Dictionary of ACMG criteria (key) and whether they are met (value)
        
        Returns:
            Tuple of (classification, criteria_met)
        """
        # Filter criteria that are met
        criteria_met = [criterion for criterion, is_met in acmg_criteria.items() if is_met]
        
        # Check for each classification, starting from most confident
        for classification, rule_sets in self.classification_rules.items():
            for rule_set in rule_sets:
                if self._check_rule_set(rule_set, criteria_met):
                    return classification, criteria_met
        
        # Default to Variant of Uncertain Significance (VUS)
        return 'uncertain_significance', criteria_met
    
    def _check_rule_set(self, rule_set: List[str], criteria_met: List[str]) -> bool:
        """
        Check if a rule set is satisfied by the given criteria.
        
        Args:
            rule_set: List of criteria patterns to check
            criteria_met: List of ACMG criteria that are met
            
        Returns:
            True if rule set is satisfied, False otherwise
        """
        for pattern in rule_set:
            # Check if at least one criterion in the pattern is met
            alternatives = pattern.split('|')
            if not any(alt in criteria_met for alt in alternatives):
                return False
        
        return True
    
    def interpret_variant(
        self, 
        variant_data: Dict, 
        gene_data: Optional[Dict] = None
    ) -> Dict:
        """
        Interpret a variant and add classification information.
        
        Args:
            variant_data: Dictionary containing variant information
            gene_data: Dictionary containing gene information (optional)
            
        Returns:
            Dictionary with additional classification information
        """
        # Extract relevant data
        result = variant_data.copy()
        
        # Determine ACMG criteria
        acmg_criteria = self._evaluate_acmg_criteria(variant_data, gene_data)
        
        # Classify variant
        classification, criteria_met = self.classify_variant(acmg_criteria)
        
        # Add classification information
        result['classification'] = self._format_classification(classification)
        result['significance_class'] = classification.replace('_', '-')
        result['acmg_criteria'] = criteria_met
        result['classification_rationale'] = self._generate_rationale(
            classification, criteria_met
        )
        
        return result
    
    def _evaluate_acmg_criteria(
        self, 
        variant_data: Dict, 
        gene_data: Optional[Dict] = None
    ) -> Dict[str, bool]:
        """
        Evaluate ACMG criteria for a variant.
        
        Args:
            variant_data: Dictionary containing variant information
            gene_data: Dictionary containing gene information (optional)
            
        Returns:
            Dictionary of ACMG criteria and whether they are met
        """
        # Initialize criteria
        criteria = {criterion: False for criterion in self.CRITERIA_WEIGHTS.keys()}
        
        # In a real implementation, this would contain complex logic to evaluate
        # each ACMG criterion based on variant properties, gene data, etc.
        # Below is a simplified example
        
        # Example: PVS1 - null variant in a gene where LOF is a known mechanism of disease
        if (
            variant_data.get('consequence', '').lower() in 
            ['frameshift', 'nonsense', 'canonical_splice'] and
            gene_data and gene_data.get('lof_mechanism', False)
        ):
            criteria['PVS1'] = True
        
        # Example: PS1 - Same amino acid change as a previously established pathogenic variant
        if variant_data.get('clinvar_same_aa', False):
            criteria['PS1'] = True
        
        # Example: PM2 - Absent from controls or at extremely low frequency
        if variant_data.get('population_af', 0) < 0.0001:
            criteria['PM2'] = True
        
        # Example: PP3 - Multiple lines of computational evidence support a deleterious effect
        if (
            variant_data.get('sift', '') == 'deleterious' and
            variant_data.get('polyphen', '') == 'probably_damaging'
        ):
            criteria['PP3'] = True
        
        # Example: BA1 - Allele frequency is > 5% in population databases
        if variant_data.get('population_af', 0) > 0.05:
            criteria['BA1'] = True
        
        # Example: BP4 - Multiple lines of computational evidence suggest no impact
        if (
            variant_data.get('sift', '') == 'tolerated' and
            variant_data.get('polyphen', '') == 'benign'
        ):
            criteria['BP4'] = True
        
        return criteria
    
    def _format_classification(self, classification: str) -> str:
        """Format classification into a human-readable string."""
        if classification == 'pathogenic':
            return 'Pathogenic'
        elif classification == 'likely_pathogenic':
            return 'Likely Pathogenic'
        elif classification == 'benign':
            return 'Benign'
        elif classification == 'likely_benign':
            return 'Likely Benign'
        else:
            return 'Uncertain Significance'
    
    def _generate_rationale(
        self, 
        classification: str, 
        criteria_met: List[str]
    ) -> str:
        """Generate rationale for variant classification."""
        if not criteria_met:
            return "Insufficient evidence for classification."
        
        criteria_descriptions = {
            'PVS1': "Null variant in a gene where loss of function is a known mechanism of disease",
            'PS1': "Same amino acid change as a previously established pathogenic variant",
            'PS2': "De novo variant in a patient with disease and no family history",
            'PS3': "Well-established functional studies show a deleterious effect",
            'PS4': "Prevalence in affected individuals significantly increased compared to controls",
            'PM1': "Located in a mutational hot spot or critical functional domain",
            'PM2': "Absent from controls or at extremely low frequency",
            'PM3': "For recessive disorders, detected in trans with a pathogenic variant",
            'PM4': "Protein length changes due to in-frame deletions/insertions or stop-loss variants",
            'PM5': "Novel missense change at an amino acid residue where a different pathogenic missense variant has been seen",
            'PM6': "Assumed de novo, but without confirmation of paternity and maternity",
            'PP1': "Co-segregation with disease in multiple affected family members",
            'PP2': "Missense variant in a gene with a low rate of benign missense variants and common pathogenic missense variants",
            'PP3': "Multiple lines of computational evidence support a deleterious effect",
            'PP4': "Patient's phenotype or family history is highly specific for a disease with a single genetic etiology",
            'PP5': "Reputable source recently reports variant as pathogenic",
            'BA1': "Allele frequency is > 5% in population databases",
            'BS1': "Allele frequency is greater than expected for disorder",
            'BS2': "Observed in healthy adult with full penetrance expected at an early age",
            'BS3': "Well-established functional studies show no deleterious effect",
            'BS4': "Lack of segregation in affected members of a family",
            'BP1': "Missense variant in a gene for which primarily truncating variants are known to cause disease",
            'BP2': "Observed in trans with a pathogenic variant for a fully penetrant dominant gene/disorder or in cis with a pathogenic variant",
            'BP3': "In-frame deletions/insertions in a repetitive region without a known function",
            'BP4': "Multiple lines of computational evidence suggest no impact on gene or gene product",
            'BP5': "Variant found in a case with an alternate molecular basis for disease",
            'BP6': "Reputable source recently reports variant as benign",
            'BP7': "A synonymous variant for which splicing prediction algorithms predict no impact"
        }
        
        # Format criteria
        criteria_text = []
        for criterion in criteria_met:
            if criterion in criteria_descriptions:
                criteria_text.append(f"{criterion}: {criteria_descriptions[criterion]}")
            else:
                criteria_text.append(criterion)
        
        # Create rationale
        if classification == 'pathogenic':
            intro = "This variant is classified as Pathogenic based on the following criteria:"
        elif classification == 'likely_pathogenic':
            intro = "This variant is classified as Likely Pathogenic based on the following criteria:"
        elif classification == 'benign':
            intro = "This variant is classified as Benign based on the following criteria:"
        elif classification == 'likely_benign':
            intro = "This variant is classified as Likely Benign based on the following criteria:"
        else:
            intro = "This variant is classified as a Variant of Uncertain Significance due to limited or conflicting evidence:"
        
        return f"{intro}\n- " + "\n- ".join(criteria_text)


def classify_variants_from_tsv(
    input_path: Union[str, Path],
    output_path: Union[str, Path],
    gene_data_path: Optional[Union[str, Path]] = None,
    rules_path: Optional[Union[str, Path]] = None
) -> None:
    """
    Classify variants from a TSV file and write results to a new TSV file.
    
    Args:
        input_path: Path to input TSV file
        output_path: Path to output TSV file
        gene_data_path: Path to gene data JSON file (optional)
        rules_path: Path to custom ACMG rules JSON file (optional)
    """
    logger.info(f"Loading variants from {input_path}")
    df = pd.read_csv(input_path, sep='\t')
    
    # Load gene data if provided
    gene_data = {}
    if gene_data_path:
        logger.info(f"Loading gene data from {gene_data_path}")
        with open(gene_data_path, 'r') as f:
            gene_data = json.load(f)
    
    # Initialize classifier
    classifier = ACMGClassifier(rules_path=rules_path)
    
    # Process variants
    logger.info("Classifying variants...")
    results = []
    
    for _, row in df.iterrows():
        # Convert row to dictionary
        variant_data = row.to_dict()
        
        # Get gene-specific data if available
        current_gene_data = gene_data.get(variant_data.get('GENE', ''), {})
        
        # Interpret variant
        result = classifier.interpret_variant(variant_data, current_gene_data)
        results.append(result)
    
    # Convert results to DataFrame
    result_df = pd.DataFrame(results)
    
    # Write to output file
    logger.info(f"Writing classified variants to {output_path}")
    result_df.to_csv(output_path, sep='\t', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ACMG variant classifier")
    parser.add_argument('--input', required=True, help='Input TSV file with variant data')
    parser.add_argument('--output', required=True, help='Output TSV file with classification')
    parser.add_argument('--gene_data', help='Gene data JSON file (optional)')
    parser.add_argument('--rules', help='Custom ACMG rules JSON file (optional)')
    
    args = parser.parse_args()
    
    try:
        classify_variants_from_tsv(
            input_path=args.input,
            output_path=args.output,
            gene_data_path=args.gene_data,
            rules_path=args.rules
        )
        logger.info("Variant classification completed successfully")
        
    except Exception as e:
        logger.error(f"Error in variant classification: {str(e)}")
        sys.exit(1) 