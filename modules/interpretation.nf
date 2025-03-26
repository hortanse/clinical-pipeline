#!/usr/bin/env nextflow

/*
========================================================================================
    INTERPRETATION MODULE
========================================================================================
    This module handles the clinical interpretation of annotated variants:
    - ACMG/AMP guideline-based variant classification
    - Phenotype-genotype correlation
    - Literature mining integration
    - Therapeutic relevance assessment
    - Reporting recommendations generation
*/

// Module parameters with defaults
params.acmg_rules = null
params.max_population_af = 0.01
params.min_coverage = 20
params.min_quality = 30
params.phenotype_data = null
params.include_clinical_notes = false

/*
 * PROCESS: ACMG/AMP classification of variants
 */
process ACMG_CLASSIFICATION {
    tag "$meta.id"
    label 'process_medium'
    container 'python:3.9'
    
    publishDir "${params.output}/interpretation", mode: 'copy'
    
    input:
    tuple val(meta), path(vcf), path(tbi), path(json)
    path acmg_rules
    path gene_data
    
    output:
    tuple val(meta), path("*.classified.json"), emit: json
    tuple val(meta), path("*.classified.vcf.gz"), emit: vcf
    tuple val(meta), path("*.classified.vcf.gz.tbi"), emit: tbi
    tuple val(meta), path("*.classification_report.html"), emit: html
    path "versions.yml", emit: versions
    
    script:
    def prefix = meta.id
    def rules_file = acmg_rules ? "--rules $acmg_rules" : ""
    def gene_data_file = gene_data ? "--gene_data $gene_data" : ""
    def max_af = params.max_population_af ? "--max_af ${params.max_population_af}" : ""
    
    """
    #!/usr/bin/env python3
    
    import json
    import gzip
    import csv
    import sys
    import os
    import re
    from pathlib import Path
    
    # This would be a more complex implementation in a real pipeline
    # that implements ACMG/AMP guidelines for variant classification
    
    # For this demonstration, we'll create a simplified version that:
    # 1. Reads the input VCF and JSON
    # 2. Applies simplified ACMG rules
    # 3. Outputs classified variants
    
    # Read the rules file
    rules = {}
    if os.path.exists("${acmg_rules}"):
        with open("${acmg_rules}") as f:
            rules = json.load(f)
    
    # Function to classify variants (simplified)
    def classify_variant(variant):
        # This would be a sophisticated algorithm in a real pipeline
        # applying ACMG/AMP rules based on variant characteristics
        
        # For demonstration, use simplified logic
        clinical_sig = variant.get("clinical_significance", "Unknown")
        consequence = variant.get("consequence", "")
        gene = variant.get("gene", "")
        af = float(variant.get("gnomad_af", "0") or "0")
        
        # Example classification logic
        if "pathogenic" in clinical_sig.lower():
            return "Pathogenic"
        elif "likely_pathogenic" in clinical_sig.lower():
            return "Likely Pathogenic"
        elif "benign" in clinical_sig.lower() and not "likely" in clinical_sig.lower():
            return "Benign"
        elif "likely_benign" in clinical_sig.lower():
            return "Likely Benign"
        elif af > ${params.max_population_af}:
            return "Likely Benign"
        elif consequence in ["frameshift", "nonsense", "splice_donor", "splice_acceptor"]:
            return "Likely Pathogenic"
        else:
            return "Uncertain Significance"
    
    # Read input JSON
    with open("${json}") as f:
        data = json.load(f)
    
    # Classify variants
    variants = data.get("variants", [])
    for variant in variants:
        variant["acmg_classification"] = classify_variant(variant)
    
    # Write classified JSON
    with open("${prefix}.classified.json", "w") as f:
        json.dump({"sample_id": "${meta.id}", "variants": variants}, f, indent=2)
    
    # Create classified VCF
    with gzip.open("${vcf}", "rt") as in_vcf, gzip.open("${prefix}.classified.vcf.gz", "wt") as out_vcf:
        for line in in_vcf:
            if line.startswith("#"):
                if line.startswith("##INFO") and "##INFO=<ID=ACMG_CLASS" not in line:
                    out_vcf.write(line)
                    if "##INFO=<ID=CLINICAL_SIG" in line:
                        out_vcf.write('##INFO=<ID=ACMG_CLASS,Number=1,Type=String,Description="ACMG/AMP classification">\n')
                else:
                    out_vcf.write(line)
            else:
                parts = line.strip().split('\\t')
                if len(parts) >= 8:
                    chrom = parts[0]
                    pos = parts[1]
                    
                    # Find matching variant in our classified list
                    matching_variant = None
                    for v in variants:
                        if v["chrom"] == chrom and v["pos"] == pos:
                            matching_variant = v
                            break
                    
                    if matching_variant:
                        # Add ACMG classification to INFO field
                        info = parts[7]
                        if not info.endswith(";"):
                            info += ";"
                        info += f"ACMG_CLASS={matching_variant['acmg_classification']}"
                        parts[7] = info
                    
                    out_vcf.write('\\t'.join(parts) + '\\n')
    
    # Index the VCF
    os.system("tabix -p vcf ${prefix}.classified.vcf.gz")
    
    # Generate an HTML report
    with open("${prefix}.classification_report.html", "w") as html_file:
        html_file.write(f"""<!DOCTYPE html>
    <html>
    <head>
        <title>Variant Classification Report - {prefix}</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 20px; }}
            h1 {{ color: #333366; }}
            table {{ border-collapse: collapse; width: 100%; }}
            th, td {{ border: 1px solid #dddddd; text-align: left; padding: 8px; }}
            tr:nth-child(even) {{ background-color: #f2f2f2; }}
            .pathogenic {{ background-color: #ffaaaa; }}
            .likely-pathogenic {{ background-color: #ffccaa; }}
            .vus {{ background-color: #ffffaa; }}
            .likely-benign {{ background-color: #aaffaa; }}
            .benign {{ background-color: #aaffcc; }}
        </style>
    </head>
    <body>
        <h1>Variant Classification Report</h1>
        <p><strong>Sample ID:</strong> {prefix}</p>
        <h2>Classified Variants</h2>
        <table>
            <tr>
                <th>Chromosome</th>
                <th>Position</th>
                <th>Reference</th>
                <th>Alternative</th>
                <th>Gene</th>
                <th>Consequence</th>
                <th>ACMG Classification</th>
            </tr>
    """)
        
        for variant in variants:
            css_class = ""
            if variant["acmg_classification"] == "Pathogenic":
                css_class = "pathogenic"
            elif variant["acmg_classification"] == "Likely Pathogenic":
                css_class = "likely-pathogenic"
            elif variant["acmg_classification"] == "Uncertain Significance":
                css_class = "vus"
            elif variant["acmg_classification"] == "Likely Benign":
                css_class = "likely-benign"
            elif variant["acmg_classification"] == "Benign":
                css_class = "benign"
            
            html_file.write(f"""
            <tr class="{css_class}">
                <td>{variant.get("chrom", "")}</td>
                <td>{variant.get("pos", "")}</td>
                <td>{variant.get("ref", "")}</td>
                <td>{variant.get("alt", "")}</td>
                <td>{variant.get("gene", "")}</td>
                <td>{variant.get("consequence", "")}</td>
                <td>{variant.get("acmg_classification", "")}</td>
            </tr>
            """)
        
        html_file.write("""
        </table>
    </body>
    </html>
    """)
    
    with open("versions.yml", "w") as versions_file:
        versions_file.write(\"""${task.process}":
        python: \$(python --version | sed 's/Python //g')
        tabix: \$(tabix --version | head -n1 | sed 's/tabix (htslib) //g')
    \""")
    """
}

/*
 * PROCESS: Phenotype correlation with genotype
 */
process PHENOTYPE_CORRELATION {
    tag "$meta.id"
    label 'process_low'
    container 'python:3.9'
    
    publishDir "${params.output}/interpretation", mode: 'copy'
    
    input:
    tuple val(meta), path(classified_json)
    path phenotype_data
    
    output:
    tuple val(meta), path("*.pheno_correlated.json"), emit: json
    tuple val(meta), path("*.phenotype_report.html"), emit: html
    path "versions.yml", emit: versions
    
    when:
    phenotype_data
    
    script:
    def prefix = meta.id
    
    """
    #!/usr/bin/env python3
    
    import json
    import csv
    import sys
    import os
    from pathlib import Path
    
    # This would be a more complex implementation in a real pipeline
    # that correlates phenotype data with genetic findings
    
    # Read classified variants
    with open("${classified_json}") as f:
        variant_data = json.load(f)
    
    # Read phenotype data
    phenotypes = []
    if os.path.exists("${phenotype_data}"):
        with open("${phenotype_data}") as f:
            if "${phenotype_data}".endswith('.json'):
                phenotype_data = json.load(f)
                phenotypes = phenotype_data.get("phenotypes", [])
            elif "${phenotype_data}".endswith('.csv'):
                reader = csv.DictReader(f)
                for row in reader:
                    if row.get("sample_id") == "${meta.id}":
                        phenotypes = row.get("phenotypes", "").split(";")
                        break
    
    # Correlate variants with phenotypes
    variants = variant_data.get("variants", [])
    for variant in variants:
        variant["phenotype_correlation"] = "Unknown"
        gene = variant.get("gene", "")
        
        # In a real implementation, this would use sophisticated algorithms
        # and databases like HPO, OMIM, etc.
        # For this demo, use simplified logic
        for phenotype in phenotypes:
            if gene in ["BRCA1", "BRCA2"] and "cancer" in phenotype.lower():
                variant["phenotype_correlation"] = "Strong"
                break
            elif gene == "CFTR" and "cystic fibrosis" in phenotype.lower():
                variant["phenotype_correlation"] = "Strong"
                break
            elif gene == "LDLR" and "cholesterol" in phenotype.lower():
                variant["phenotype_correlation"] = "Strong"
                break
    
    # Write correlated JSON
    with open("${prefix}.pheno_correlated.json", "w") as f:
        json.dump({"sample_id": "${meta.id}", "variants": variants, "phenotypes": phenotypes}, f, indent=2)
    
    # Generate an HTML report
    with open("${prefix}.phenotype_report.html", "w") as html_file:
        html_file.write(f"""<!DOCTYPE html>
    <html>
    <head>
        <title>Phenotype Correlation Report - {prefix}</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 20px; }}
            h1 {{ color: #333366; }}
            table {{ border-collapse: collapse; width: 100%; }}
            th, td {{ border: 1px solid #dddddd; text-align: left; padding: 8px; }}
            tr:nth-child(even) {{ background-color: #f2f2f2; }}
            .strong {{ background-color: #aaddff; }}
        </style>
    </head>
    <body>
        <h1>Phenotype Correlation Report</h1>
        <p><strong>Sample ID:</strong> {prefix}</p>
        
        <h2>Reported Phenotypes</h2>
        <ul>
    """)
        
        for phenotype in phenotypes:
            html_file.write(f"<li>{phenotype}</li>\\n")
        
        html_file.write("""
        </ul>
        
        <h2>Genotype-Phenotype Correlations</h2>
        <table>
            <tr>
                <th>Gene</th>
                <th>Variant</th>
                <th>Classification</th>
                <th>Phenotype Correlation</th>
            </tr>
        """)
        
        for variant in variants:
            css_class = "strong" if variant.get("phenotype_correlation") == "Strong" else ""
            html_file.write(f"""
            <tr class="{css_class}">
                <td>{variant.get("gene", "")}</td>
                <td>{variant.get("chrom", "")}:{variant.get("pos", "")} {variant.get("ref", "")}&gt;{variant.get("alt", "")}</td>
                <td>{variant.get("acmg_classification", "")}</td>
                <td>{variant.get("phenotype_correlation", "")}</td>
            </tr>
            """)
        
        html_file.write("""
        </table>
    </body>
    </html>
    """)
    
    with open("versions.yml", "w") as versions_file:
        versions_file.write(\"""${task.process}":
        python: \$(python --version | sed 's/Python //g')
    \""")
    """
}

/*
 * PROCESS: Clinical relevance assessment
 */
process CLINICAL_RELEVANCE {
    tag "$meta.id"
    label 'process_low'
    container 'python:3.9'
    
    publishDir "${params.output}/interpretation", mode: 'copy'
    
    input:
    tuple val(meta), path(correlated_json)
    
    output:
    tuple val(meta), path("*.final_interpretation.json"), emit: json
    tuple val(meta), path("*.clinical_report.html"), emit: html
    path "versions.yml", emit: versions
    
    script:
    def prefix = meta.id
    
    """
    #!/usr/bin/env python3
    
    import json
    import datetime
    
    # Read correlated variants
    with open("${correlated_json}") as f:
        data = json.load(f)
    
    # Add clinical relevance assessment
    variants = data.get("variants", [])
    phenotypes = data.get("phenotypes", [])
    
    # Find clinically significant variants
    significant_variants = []
    for variant in variants:
        classification = variant.get("acmg_classification", "")
        correlation = variant.get("phenotype_correlation", "")
        
        if classification in ["Pathogenic", "Likely Pathogenic"]:
            variant["clinical_relevance"] = "High"
            significant_variants.append(variant)
        elif classification == "Uncertain Significance" and correlation == "Strong":
            variant["clinical_relevance"] = "Moderate"
            significant_variants.append(variant)
        else:
            variant["clinical_relevance"] = "Low"
    
    # Generate clinical recommendation
    recommendation = ""
    if significant_variants:
        recommendation = "Based on the genomic analysis, the following clinical actions are recommended:\\n"
        for variant in significant_variants:
            gene = variant.get("gene", "")
            if gene in ["BRCA1", "BRCA2"]:
                recommendation += f"- Referral to cancer genetics due to {gene} {variant.get('acmg_classification', '')} variant\\n"
            elif gene == "LDLR":
                recommendation += "- Lipid profile monitoring and consideration of statin therapy\\n"
            elif gene == "CFTR":
                recommendation += "- Pulmonary function testing and sweat chloride test\\n"
            else:
                recommendation += f"- Clinical correlation for {gene} variant\\n"
    else:
        recommendation = "No clinically significant findings requiring immediate action were identified."
    
    # Add recommendation to data
    data["clinical_recommendation"] = recommendation
    data["report_date"] = datetime.datetime.now().strftime("%Y-%m-%d")
    
    # Write final interpretation JSON
    with open("${prefix}.final_interpretation.json", "w") as f:
        json.dump(data, f, indent=2)
    
    # Generate an HTML report
    with open("${prefix}.clinical_report.html", "w") as html_file:
        html_file.write(f"""<!DOCTYPE html>
    <html>
    <head>
        <title>Clinical Relevance Report - {prefix}</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 20px; }}
            h1 {{ color: #333366; }}
            h2 {{ color: #336699; }}
            table {{ border-collapse: collapse; width: 100%; margin-bottom: 20px; }}
            th, td {{ border: 1px solid #dddddd; text-align: left; padding: 8px; }}
            tr:nth-child(even) {{ background-color: #f2f2f2; }}
            .high {{ background-color: #ffaaaa; }}
            .moderate {{ background-color: #ffffaa; }}
            .recommendation {{ background-color: #eeeeff; padding: 10px; border-left: 5px solid #336699; }}
            .date {{ color: #666666; font-style: italic; }}
        </style>
    </head>
    <body>
        <h1>Clinical Relevance Assessment</h1>
        <p><strong>Sample ID:</strong> {prefix}</p>
        <p class="date"><strong>Report Date:</strong> {data["report_date"]}</p>
        
        <h2>Clinically Significant Findings</h2>
        <table>
            <tr>
                <th>Gene</th>
                <th>Variant</th>
                <th>Classification</th>
                <th>Phenotype Correlation</th>
                <th>Clinical Relevance</th>
            </tr>
    """)
        
        for variant in variants:
            relevance = variant.get("clinical_relevance", "")
            if relevance in ["High", "Moderate"]:
                css_class = relevance.lower()
                html_file.write(f"""
                <tr class="{css_class}">
                    <td>{variant.get("gene", "")}</td>
                    <td>{variant.get("chrom", "")}:{variant.get("pos", "")} {variant.get("ref", "")}&gt;{variant.get("alt", "")}</td>
                    <td>{variant.get("acmg_classification", "")}</td>
                    <td>{variant.get("phenotype_correlation", "")}</td>
                    <td>{relevance}</td>
                </tr>
                """)
        
        html_file.write("""
        </table>
        
        <h2>Clinical Recommendations</h2>
        <div class="recommendation">
        """)
        
        for line in recommendation.split('\\n'):
            if line:
                html_file.write(f"<p>{line}</p>\\n")
        
        html_file.write("""
        </div>
        
        <h2>Methodology</h2>
        <p>Variants were classified according to ACMG/AMP guidelines and correlated with reported phenotypes. 
        Clinical relevance was assessed based on variant classification, phenotype correlation, and known disease associations.</p>
    </body>
    </html>
    """)
    
    with open("versions.yml", "w") as versions_file:
        versions_file.write(\"""${task.process}":
        python: \$(python --version | sed 's/Python //g')
    \""")
    """
}

// Main workflow for interpretation
workflow INTERPRETATION {
    take:
    vcf_tbi_json         // [meta, vcf, tbi, json]
    acmg_rules           // ACMG classification rules
    gene_data            // Gene-disease associations and other gene data
    phenotype_data       // Optional: phenotype data for correlation
    
    main:
    ch_versions = Channel.empty()
    
    // Classify variants using ACMG/AMP guidelines
    ACMG_CLASSIFICATION(
        vcf_tbi_json,
        acmg_rules,
        gene_data
    )
    ch_versions = ch_versions.mix(ACMG_CLASSIFICATION.out.versions)
    
    // Correlate with phenotypes if available
    if (phenotype_data) {
        PHENOTYPE_CORRELATION(
            ACMG_CLASSIFICATION.out.json,
            phenotype_data
        )
        ch_versions = ch_versions.mix(PHENOTYPE_CORRELATION.out.versions)
        
        CLINICAL_RELEVANCE(
            PHENOTYPE_CORRELATION.out.json
        )
        
        ch_final_json = CLINICAL_RELEVANCE.out.json
        ch_clinical_html = CLINICAL_RELEVANCE.out.html
        ch_versions = ch_versions.mix(CLINICAL_RELEVANCE.out.versions)
    } else {
        // Skip phenotype correlation if no phenotype data
        CLINICAL_RELEVANCE(
            ACMG_CLASSIFICATION.out.json
        )
        
        ch_final_json = CLINICAL_RELEVANCE.out.json
        ch_clinical_html = CLINICAL_RELEVANCE.out.html
        ch_versions = ch_versions.mix(CLINICAL_RELEVANCE.out.versions)
    }
    
    emit:
    classified_vcf = ACMG_CLASSIFICATION.out.vcf
    classified_tbi = ACMG_CLASSIFICATION.out.tbi
    classified_json = ACMG_CLASSIFICATION.out.json
    final_json = ch_final_json
    classification_report = ACMG_CLASSIFICATION.out.html
    clinical_report = ch_clinical_html
    versions = ch_versions
} 