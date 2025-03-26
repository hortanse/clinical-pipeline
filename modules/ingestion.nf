/*
========================================================================================
    Sample Ingestion Module
========================================================================================
    This module handles sample accession and FASTQ retrieval from LIMS
*/

process FETCH_SAMPLES {
    tag "$samplesheet"
    label 'process_low'
    
    input:
    path samplesheet
    
    output:
    path "samples.json", emit: samples_json
    path "fastq/*{1,2}.fastq.gz", emit: reads
    
    script:
...(about 7 lines omitted)...
process VALIDATE_SAMPLES {
    tag "$samples_json"
    label 'process_low'
    
    input:
    path samples_json
    
    output:
    path "validated_samples.json", emit: validated_samples
    path "validation_report.html", emit: validation_report
    
    script:
    """
    python ${projectDir}/lib/lims/validate_samples.py \\
        --samples ${samples_json} \\
        --output validated_samples.json \\
        --report validation_report.html
    """
} 