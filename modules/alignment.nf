/*
========================================================================================
    Alignment Module
========================================================================================
    This module handles read alignment to reference genome and post-alignment processing
*/

process ALIGN_READS {
    tag "${meta.id}"
    label 'process_high'
    
    input:
    tuple val(meta), path(reads)
    path genome
    
    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.bai"), emit: bai
    path "versions.yml", emit: versions
    
    script:
    def prefix = "${meta.id}"
    def read_group = "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tLB:${meta.library ?: 'unknown'}\\tPL:ILLUMINA"
    
    """
    # Align reads to reference genome with BWA MEM
    bwa mem \\
        -t $task.cpus \\
        -R \"${read_group}\" \\
        ${genome}/BWAIndex/genome.fa \\
        ${reads[0]} ${reads[1]} \\
        | samtools view -@ $task.cpus -Sbh - \\
        | samtools sort -@ $task.cpus -o ${prefix}.sorted.bam -
    
    # Index BAM file
    samtools index ${prefix}.sorted.bam ${prefix}.sorted.bai
    
    # Rename files to ensure consistent naming convention
    mv ${prefix}.sorted.bam ${prefix}.bam
    mv ${prefix}.sorted.bai ${prefix}.bai
    
    # Record software versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | grep -e '^Version' | sed 's/Version: //g')
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //g')
    END_VERSIONS
    """
}

process MARK_DUPLICATES {
    tag "${meta.id}"
    label 'process_medium'
    
    input:
    tuple val(meta), path(bam)
    
    output:
    tuple val(meta), path("*.markdup.bam"), emit: bam
    tuple val(meta), path("*.markdup.bai"), emit: bai
    path "versions.yml", emit: versions
    
    script:
    def prefix = "${meta.id}"
    def avail_mem = task.memory ? "-Xmx${task.memory.toGiga()}g" : ""
    def known_sites_command = known_sites.collect{ "--known-sites $it" }.join(' ')
    
    """
    # Generate recalibration table
    gatk --java-options "${avail_mem}" BaseRecalibrator \\
        -I ${bam} \\
        -R ${genome}/genome.fa \\
        ${known_sites_command} \\
        -O ${prefix}.recal.table
    
    # Apply base recalibration
    gatk --java-options "${avail_mem}" ApplyBQSR \\
        -I ${bam} \\
        -R ${genome}/genome.fa \\
        --bqsr-recal-file ${prefix}.recal.table \\
        -O ${prefix}.recal.bam
    
    # Record software versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version | grep -e "^GATK" | sed 's/GATK //g')
    END_VERSIONS
    """
}

process COLLECT_ALIGNMENT_METRICS {
    tag "${meta.id}"
    label 'process_low'
    
    input:
    tuple val(meta), path(bam), path(bai)
    path genome
    
    output:
    tuple val(meta), path("*.alignment_metrics.txt"), emit: metrics
    path "versions.yml", emit: versions
    
    script:
    def prefix = "${meta.id}"
    def avail_mem = task.memory ? "-Xmx${task.memory.toGiga()}g" : ""
    
    """
    # Collect alignment metrics
    gatk --java-options "${avail_mem}" CollectAlignmentSummaryMetrics \\
        -R ${genome}/genome.fa \\
        -I ${bam} \\
        -O ${prefix}.alignment_metrics.txt
    
    # Record software versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version | grep -e "^GATK" | sed 's/GATK //g')
    END_VERSIONS
    """
} 