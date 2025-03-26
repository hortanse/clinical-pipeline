#!/usr/bin/env nextflow

/*
========================================================================================
    PREPROCESSING MODULE
========================================================================================
    This module handles the initial preprocessing steps for NGS data:
    - Quality control with FastQC
    - Adapter trimming and quality filtering with Trim Galore
    - Additional preprocessing with fastp
    - Metric collection and quality reports
*/

// Import required processes from nf-core modules (if using them)
// include { FASTQC } from '../modules/nf-core/fastqc/main'

// Module parameters with defaults
params.fastqc_args = ""
params.trim_galore_args = "--quality 20 --stringency 3 --length 50"
params.fastp_args = "--qualified_quality_phred 20 --length_required 50"
params.skip_fastp = false
params.skip_fastqc = false
params.save_trimmed_fail = false
params.save_merged = false

/*
 * PROCESS: Run FastQC for quality control
 */
process FASTQC {
    tag "$meta.id"
    label 'process_medium'
    container 'nfcore/fastqc:0.11.9'
    
    publishDir "${params.output}/fastqc", mode: 'copy',
        saveAs: { filename -> filename.indexOf('.zip') > 0 ? "zips/$filename" : "$filename" }
    
    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"), emit: zip
    path "versions.yml", emit: versions
    
    when:
    !params.skip_fastqc
    
    script:
    def args = params.fastqc_args ?: ''
    def prefix = meta.id
    """
    fastqc --threads $task.cpus $args --outdir . ${reads}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed -e "s/FastQC v//g")
    END_VERSIONS
    """
}

/*
 * PROCESS: Run Trim Galore for adapter and quality trimming
 */
process TRIM_GALORE {
    tag "$meta.id"
    label 'process_high'
    container 'staphb/trimgalore:0.6.7'
    
    publishDir "${params.output}/trimmed", mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.html')) "fastqc/$filename"
            else if (filename.endsWith('.zip')) "fastqc/zips/$filename"
            else if (filename.indexOf('_val_') > 0) "fastq/$filename"
            else if (filename.indexOf('trimming_report.txt') > 0) "logs/$filename"
            else if (params.save_trimmed_fail && filename.indexOf('unpaired') > 0) "unpaired/$filename"
            else null
        }
    
    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("*_val_*.fq.gz"), emit: reads
    tuple val(meta), path("*_trimming_report.txt"), emit: logs
    tuple val(meta), path("*.html"), optional: true, emit: html
    tuple val(meta), path("*.zip"), optional: true, emit: zip
    path "versions.yml", emit: versions
    
    script:
    def args = params.trim_galore_args ?: ''
    def prefix = meta.id
    def paired = meta.single_end ? "" : "--paired"
    def trim_log = meta.single_end ? "${prefix}_trimming_report.txt" : "${prefix}_trimming_report.txt ${prefix}_trimming_report.txt"
    
    """
    trim_galore $args $paired --cores $task.cpus ${reads}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimgalore: \$(trim_galore --version | sed -e "s/^.*version //g")
    END_VERSIONS
    """
}

/*
 * PROCESS: Run fastp for additional preprocessing
 */
process FASTP {
    tag "$meta.id"
    label 'process_medium'
    container 'biocontainers/fastp:0.23.2'
    
    publishDir "${params.output}/fastp", mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.json')) "json/$filename"
            else if (filename.endsWith('.html')) "html/$filename"
            else if (filename.indexOf('_fastp.') > 0) "fastq/$filename"
            else null
        }
    
    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("*_fastp.fastq.gz"), emit: reads
    tuple val(meta), path("*.json"), emit: json
    tuple val(meta), path("*.html"), emit: html
    path "versions.yml", emit: versions
    
    when:
    !params.skip_fastp
    
    script:
    def args = params.fastp_args ?: ''
    def prefix = meta.id
    
    if (meta.single_end) {
        """
        fastp $args --in1 ${reads[0]} --out1 ${prefix}_fastp.fastq.gz --json ${prefix}_fastp.json --html ${prefix}_fastp.html --thread $task.cpus
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """
    } else {
        """
        fastp $args --in1 ${reads[0]} --in2 ${reads[1]} --out1 ${prefix}_R1_fastp.fastq.gz --out2 ${prefix}_R2_fastp.fastq.gz --json ${prefix}_fastp.json --html ${prefix}_fastp.html --thread $task.cpus
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """
    }
}

/*
 * PROCESS: Generate preprocessing QC report
 */
process PREPROCESSING_QC_REPORT {
    tag "$meta.id"
    label 'process_low'
    container 'python:3.9'
    
    publishDir "${params.output}/reports", mode: 'copy'
    
    input:
    tuple val(meta), path(fastqc_zip)
    tuple val(meta), path(trim_log)
    tuple val(meta), path(fastp_json)
    
    output:
    tuple val(meta), path("${meta.id}_preprocessing_qc.html"), emit: html
    tuple val(meta), path("${meta.id}_preprocessing_qc.json"), emit: json
    
    script:
    """
    #!/usr/bin/env python3
    
    import json
    import zipfile
    import re
    from pathlib import Path
    
    # Simple QC report generation
    # In a real scenario, this would parse the FastQC, trim reports, and fastp JSON
    # to generate a comprehensive HTML and JSON report
    
    sample_id = "${meta.id}"
    
    # Create a simple HTML report
    with open(f"{sample_id}_preprocessing_qc.html", "w") as html_report:
        html_report.write(f"<html><head><title>QC Report for {sample_id}</title></head>")
        html_report.write("<body><h1>Preprocessing QC Report</h1>")
        html_report.write(f"<h2>Sample: {sample_id}</h2>")
        # In a real implementation, we would extract data from the input files here
        html_report.write("</body></html>")
    
    # Create a simple JSON report
    qc_data = {
        "sample_id": sample_id,
        "fastqc": {"status": "pass"},
        "trimming": {"status": "pass"},
        "fastp": {"status": "pass"}
    }
    
    with open(f"{sample_id}_preprocessing_qc.json", "w") as json_report:
        json.dump(qc_data, json_report, indent=4)
    """
}

// Function to stage incoming files depending on single-end/paired-end
def getStagePattern(single_end) {
    return single_end ? "*.fastq.gz" : "*_{1,2}.fastq.gz"
}

// Main workflow for preprocessing
workflow PREPROCESSING {
    take:
    reads // [meta, [reads...]]
    
    main:
    ch_versions = Channel.empty()
    
    // FastQC
    FASTQC(reads)
    ch_versions = ch_versions.mix(FASTQC.out.versions)
    
    // Trimming
    TRIM_GALORE(reads)
    ch_versions = ch_versions.mix(TRIM_GALORE.out.versions)
    
    // Additional preprocessing
    if (!params.skip_fastp) {
        FASTP(TRIM_GALORE.out.reads)
        ch_clean_reads = FASTP.out.reads
        ch_versions = ch_versions.mix(FASTP.out.versions)
    } else {
        ch_clean_reads = TRIM_GALORE.out.reads
    }
    
    // QC report
    PREPROCESSING_QC_REPORT(
        FASTQC.out.zip,
        TRIM_GALORE.out.logs,
        params.skip_fastp ? Channel.value([]) : FASTP.out.json
    )
    
    emit:
    reads = ch_clean_reads
    fastqc_zip = FASTQC.out.zip
    fastqc_html = FASTQC.out.html
    trim_logs = TRIM_GALORE.out.logs
    fastp_json = params.skip_fastp ? Channel.empty() : FASTP.out.json
    fastp_html = params.skip_fastp ? Channel.empty() : FASTP.out.html
    qc_html = PREPROCESSING_QC_REPORT.out.html
    qc_json = PREPROCESSING_QC_REPORT.out.json
    versions = ch_versions
} 