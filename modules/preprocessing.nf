/*
========================================================================================
    Quality Control and Preprocessing Module
========================================================================================
    This module handles initial QC and read preprocessing steps
*/

process FASTQC {
    tag "${meta.id}"
    label 'process_medium'
    
    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("*_fastqc.html"), emit: html
    tuple val(meta), path("*_fastqc.zip"), emit: zip
    path "versions.yml", emit: versions
    
    script:
...(about 58 lines omitted)...
    """
    fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${prefix}_fastp_1.fastq.gz \\
        --out2 ${prefix}_fastp_2.fastq.gz \\
        --json ${prefix}_fastp.json \\
        --html ${prefix}_fastp.html \\
        --thread $task.cpus \\
        --detect_adapter_for_pe \\
        --qualified_quality_phred 20 \\
        --length_required 50
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/^.*fastp //g")
    END_VERSIONS
    """
} 