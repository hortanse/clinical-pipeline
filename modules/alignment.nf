#!/usr/bin/env nextflow

/*
========================================================================================
    ALIGNMENT MODULE
========================================================================================
    This module handles the alignment of processed reads to the reference genome:
    - BWA-MEM alignment
    - SAM to BAM conversion
    - Sorting and indexing
    - Duplicate marking
    - Quality metric collection
    - BAM file merge for multiple libraries
*/

// Module parameters with defaults
params.bwa_args = ""
params.markdup_java_options = "-Xmx8g"
params.save_unaligned = false
params.save_merged_bam = true
params.save_unmapped = false
params.sort_bam = true

/*
 * PROCESS: Align reads with BWA-MEM
 */
process BWA_MEM {
    tag "$meta.id"
    label 'process_high'
    container 'biocontainers/bwa:v0.7.17_cv1'
    
    input:
    tuple val(meta), path(reads)
    path(reference)
    path(reference_index) // BWA index files
    
    output:
    tuple val(meta), path("*.sam"), emit: sam
    path "versions.yml", emit: versions
    
    script:
    def args = params.bwa_args ?: ''
    def prefix = meta.id
    def read_list = reads.collect {it.toString()}.join(' ')
    
    """
    bwa mem \\
        -t $task.cpus \\
        $args \\
        -R '@RG\\tID:${prefix}\\tSM:${prefix}\\tPL:ILLUMINA' \\
        $reference \\
        $read_list > ${prefix}.sam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | grep 'Version' | sed 's/Version: //')
    END_VERSIONS
    """
}

/*
 * PROCESS: Convert SAM to BAM
 */
process SAM_TO_BAM {
    tag "$meta.id"
    label 'process_medium'
    container 'biocontainers/samtools:1.15.1'
    
    input:
    tuple val(meta), path(sam)
    
    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml", emit: versions
    
    script:
    def prefix = meta.id
    
    """
    samtools view -@ $task.cpus -bS $sam > ${prefix}.bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | grep samtools | sed 's/samtools //')
    END_VERSIONS
    """
}

/*
 * PROCESS: Sort BAM file
 */
process SORT_BAM {
    tag "$meta.id"
    label 'process_medium'
    container 'biocontainers/samtools:1.15.1'
    
    publishDir "${params.output}/aligned", mode: 'copy',
        saveAs: { filename -> params.save_merged_bam && meta.id == meta.sample_id ? filename : null }
    
    input:
    tuple val(meta), path(bam)
    
    output:
    tuple val(meta), path("*.sorted.bam"), emit: bam
    tuple val(meta), path("*.sorted.bam.bai"), emit: bai
    path "versions.yml", emit: versions
    
    when:
    params.sort_bam
    
    script:
    def prefix = meta.id
    
    """
    samtools sort -@ $task.cpus -o ${prefix}.sorted.bam $bam
    samtools index -@ $task.cpus ${prefix}.sorted.bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | grep samtools | sed 's/samtools //')
    END_VERSIONS
    """
}

/*
 * PROCESS: Mark duplicates
 */
process MARK_DUPLICATES {
    tag "$meta.id"
    label 'process_medium'
    container 'broadinstitute/gatk:4.3.0.0'
    
    publishDir "${params.output}/aligned", mode: 'copy',
        saveAs: { filename -> 
            if (filename.endsWith('.metrics.txt')) "metrics/$filename"
            else if (params.save_merged_bam && meta.id == meta.sample_id) filename
            else null
        }
    
    input:
    tuple val(meta), path(bam)
    
    output:
    tuple val(meta), path("*.markdup.bam"), emit: bam
    tuple val(meta), path("*.markdup.bam.bai"), emit: bai
    tuple val(meta), path("*.metrics.txt"), emit: metrics
    path "versions.yml", emit: versions
    
    script:
    def prefix = meta.id
    def java_opts = params.markdup_java_options ?: ''
    
    """
    gatk --java-options "${java_opts}" MarkDuplicates \\
        -I $bam \\
        -O ${prefix}.markdup.bam \\
        -M ${prefix}.metrics.txt \\
        --CREATE_INDEX true
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version | grep GATK | sed 's/GATK v//')
    END_VERSIONS
    """
}

/*
 * PROCESS: Merge BAM files for the same sample
 */
process MERGE_BAM {
    tag "$meta.id"
    label 'process_medium'
    container 'biocontainers/samtools:1.15.1'
    
    publishDir "${params.output}/aligned", mode: 'copy'
    
    input:
    tuple val(meta), path(bams)
    
    output:
    tuple val(meta), path("*_merged.bam"), emit: bam
    tuple val(meta), path("*_merged.bam.bai"), emit: bai
    path "versions.yml", emit: versions
    
    script:
    def prefix = meta.id
    def bam_files = bams instanceof List ? bams.join(' ') : bams
    
    """
    samtools merge -@ $task.cpus ${prefix}_merged.bam $bam_files
    samtools index -@ $task.cpus ${prefix}_merged.bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | grep samtools | sed 's/samtools //')
    END_VERSIONS
    """
}

/*
 * PROCESS: Collect alignment metrics
 */
process ALIGNMENT_METRICS {
    tag "$meta.id"
    label 'process_low'
    container 'broadinstitute/gatk:4.3.0.0'
    
    publishDir "${params.output}/aligned/metrics", mode: 'copy'
    
    input:
    tuple val(meta), path(bam), path(bai)
    path(reference)
    path(reference_dict)
    
    output:
    tuple val(meta), path("*.alignment_metrics.txt"), emit: metrics
    tuple val(meta), path("*.insert_size_metrics.txt"), emit: insert_metrics
    tuple val(meta), path("*.insert_size_histogram.pdf"), emit: insert_histogram
    path "versions.yml", emit: versions
    
    script:
    def prefix = meta.id
    
    """
    gatk CollectAlignmentSummaryMetrics \\
        -R $reference \\
        -I $bam \\
        -O ${prefix}.alignment_metrics.txt
    
    gatk CollectInsertSizeMetrics \\
        -I $bam \\
        -O ${prefix}.insert_size_metrics.txt \\
        -H ${prefix}.insert_size_histogram.pdf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version | grep GATK | sed 's/GATK v//')
    END_VERSIONS
    """
}

// Main workflow for alignment
workflow ALIGNMENT {
    take:
    reads              // [meta, [reads...]]
    reference          // reference fasta
    reference_index    // bwa index files
    reference_dict     // sequence dictionary
    
    main:
    ch_versions = Channel.empty()
    
    // Align reads with BWA-MEM
    BWA_MEM(reads, reference, reference_index)
    ch_versions = ch_versions.mix(BWA_MEM.out.versions)
    
    // Convert SAM to BAM
    SAM_TO_BAM(BWA_MEM.out.sam)
    ch_versions = ch_versions.mix(SAM_TO_BAM.out.versions)
    
    // Sort BAM file
    SORT_BAM(SAM_TO_BAM.out.bam)
    ch_versions = ch_versions.mix(SORT_BAM.out.versions)
    
    // Mark duplicates
    MARK_DUPLICATES(SORT_BAM.out.bam)
    ch_versions = ch_versions.mix(MARK_DUPLICATES.out.versions)
    
    // Group BAMs by sample ID for potential merging
    ch_bam_to_merge = MARK_DUPLICATES.out.bam
        .map { meta, bam -> 
            def new_meta = [id: meta.sample_id, sample_id: meta.sample_id]
            [new_meta, bam] 
        }
        .groupTuple()
    
    // Merge BAMs if there are multiple per sample
    MERGE_BAM(ch_bam_to_merge)
    ch_versions = ch_versions.mix(MERGE_BAM.out.versions)
    
    // Collect alignment metrics
    ALIGNMENT_METRICS(
        MERGE_BAM.out.bam.join(MERGE_BAM.out.bai, by: [0]), 
        reference,
        reference_dict
    )
    ch_versions = ch_versions.mix(ALIGNMENT_METRICS.out.versions)
    
    emit:
    bam = MERGE_BAM.out.bam
    bai = MERGE_BAM.out.bai
    metrics = ALIGNMENT_METRICS.out.metrics
    insert_metrics = ALIGNMENT_METRICS.out.insert_metrics
    versions = ch_versions
} 