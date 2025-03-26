#!/usr/bin/env nextflow

/*
========================================================================================
    VARIANT CALLING MODULE
========================================================================================
    This module handles variant calling from aligned reads:
    - GATK HaplotypeCaller for germline variant calling
    - GATK Mutect2 for somatic variant calling
    - GATK VariantFiltration for basic filtering
    - GATK VQSR for variant quality score recalibration
    - Optional: structural variant calling
*/

// Module parameters with defaults
params.variant_caller = "haplotypecaller"  // Options: haplotypecaller, mutect2
params.haplotypecaller_args = "--standard-min-confidence-threshold-for-calling 20.0"
params.mutect2_args = "--f1r2-tar-gz f1r2.tar.gz"
params.apply_vqsr = true
params.vqsr_snp_args = "--trust-all-polymorphic"
params.vqsr_indel_args = "--trust-all-polymorphic"
params.save_filtered_variants = true
params.save_gvcf = false

/*
 * PROCESS: GATK HaplotypeCaller for germline variant calling
 */
process HAPLOTYPECALLER {
    tag "$meta.id"
    label 'process_high'
    container 'broadinstitute/gatk:4.3.0.0'
    
    publishDir "${params.output}/variants/raw", mode: 'copy',
        saveAs: { filename -> 
            if (filename.endsWith('.g.vcf.gz')) "gvcf/$filename"
            else if (filename.endsWith('.vcf.gz')) "$filename"
            else null
        }
    
    input:
    tuple val(meta), path(bam), path(bai)
    path(reference)
    path(reference_dict)
    path(reference_fai)
    path(regions) // Optional: interval list for targeted sequencing
    
    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    tuple val(meta), path("*.g.vcf.gz"), optional: true, emit: gvcf
    tuple val(meta), path("*.g.vcf.gz.tbi"), optional: true, emit: gvcf_tbi
    path "versions.yml", emit: versions
    
    when:
    params.variant_caller == "haplotypecaller"
    
    script:
    def args = params.haplotypecaller_args ?: ''
    def prefix = meta.id
    def interval_command = regions ? "--intervals $regions" : ""
    def gvcf_command = params.save_gvcf ? "-ERC GVCF" : ""
    
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" HaplotypeCaller \\
        -R $reference \\
        -I $bam \\
        $interval_command \\
        $gvcf_command \\
        -O ${prefix}.vcf.gz \\
        $args
    
    if [ "$gvcf_command" != "" ]; then
        gatk --java-options "-Xmx${task.memory.toGiga()}g" HaplotypeCaller \\
            -R $reference \\
            -I $bam \\
            $interval_command \\
            -ERC GVCF \\
            -O ${prefix}.g.vcf.gz \\
            $args
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version | grep GATK | sed 's/GATK v//')
    END_VERSIONS
    """
}

/*
 * PROCESS: GATK Mutect2 for somatic variant calling
 */
process MUTECT2 {
    tag "$meta.id"
    label 'process_high'
    container 'broadinstitute/gatk:4.3.0.0'
    
    publishDir "${params.output}/variants/raw", mode: 'copy'
    
    input:
    tuple val(meta), path(bam), path(bai)
    path(reference)
    path(reference_dict)
    path(reference_fai)
    path(regions) // Optional: interval list for targeted sequencing
    tuple val(normal_meta), path(normal_bam), path(normal_bai), optional: true // Optional: matched normal
    
    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    tuple val(meta), path("*.f1r2.tar.gz"), emit: f1r2
    path "versions.yml", emit: versions
    
    when:
    params.variant_caller == "mutect2"
    
    script:
    def args = params.mutect2_args ?: ''
    def prefix = meta.id
    def interval_command = regions ? "--intervals $regions" : ""
    def normal_command = normal_bam ? "--normal-sample ${normal_meta.id} --input ${normal_bam}" : ""
    
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" Mutect2 \\
        -R $reference \\
        -I $bam \\
        -tumor ${meta.id} \\
        $normal_command \\
        $interval_command \\
        -O ${prefix}.vcf.gz \\
        --f1r2-tar-gz ${prefix}.f1r2.tar.gz \\
        $args
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version | grep GATK | sed 's/GATK v//')
    END_VERSIONS
    """
}

/*
 * PROCESS: GATK VariantFiltration for simple filtering
 */
process VARIANT_FILTRATION {
    tag "$meta.id"
    label 'process_medium'
    container 'broadinstitute/gatk:4.3.0.0'
    
    publishDir "${params.output}/variants/filtered", mode: 'copy',
        saveAs: { filename -> 
            if (params.save_filtered_variants) filename
            else null
        }
    
    input:
    tuple val(meta), path(vcf), path(tbi)
    path(reference)
    path(reference_dict)
    path(reference_fai)
    
    output:
    tuple val(meta), path("*.filtered.vcf.gz"), emit: vcf
    tuple val(meta), path("*.filtered.vcf.gz.tbi"), emit: tbi
    path "versions.yml", emit: versions
    
    when:
    !params.apply_vqsr
    
    script:
    def prefix = meta.id
    
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" VariantFiltration \\
        -R $reference \\
        -V $vcf \\
        --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \\
        --filter-name "basic_snp_filter" \\
        -O ${prefix}.filtered.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version | grep GATK | sed 's/GATK v//')
    END_VERSIONS
    """
}

/*
 * PROCESS: GATK VariantRecalibrator for SNPs
 */
process VARIANT_RECALIBRATOR_SNP {
    tag "$meta.id"
    label 'process_high'
    container 'broadinstitute/gatk:4.3.0.0'
    
    input:
    tuple val(meta), path(vcf), path(tbi)
    path(reference)
    path(reference_dict)
    path(reference_fai)
    path(hapmap)
    path(hapmap_idx)
    path(omni)
    path(omni_idx)
    path(thousandG)
    path(thousandG_idx)
    path(dbsnp)
    path(dbsnp_idx)
    
    output:
    tuple val(meta), path("*.snp.recal"), emit: recal
    tuple val(meta), path("*.snp.tranches"), emit: tranches
    tuple val(meta), path("*.snp.plots.R"), emit: plots
    path "versions.yml", emit: versions
    
    when:
    params.apply_vqsr
    
    script:
    def prefix = meta.id
    def args = params.vqsr_snp_args ?: ''
    
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" VariantRecalibrator \\
        -R $reference \\
        -V $vcf \\
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \\
        --resource:omni,known=false,training=true,truth=false,prior=12.0 $omni \\
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 $thousandG \\
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \\
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \\
        -mode SNP \\
        -O ${prefix}.snp.recal \\
        --tranches-file ${prefix}.snp.tranches \\
        --rscript-file ${prefix}.snp.plots.R \\
        $args
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version | grep GATK | sed 's/GATK v//')
    END_VERSIONS
    """
}

/*
 * PROCESS: GATK VariantRecalibrator for INDELs
 */
process VARIANT_RECALIBRATOR_INDEL {
    tag "$meta.id"
    label 'process_high'
    container 'broadinstitute/gatk:4.3.0.0'
    
    input:
    tuple val(meta), path(vcf), path(tbi)
    path(reference)
    path(reference_dict)
    path(reference_fai)
    path(mills)
    path(mills_idx)
    path(dbsnp)
    path(dbsnp_idx)
    
    output:
    tuple val(meta), path("*.indel.recal"), emit: recal
    tuple val(meta), path("*.indel.tranches"), emit: tranches
    tuple val(meta), path("*.indel.plots.R"), emit: plots
    path "versions.yml", emit: versions
    
    when:
    params.apply_vqsr
    
    script:
    def prefix = meta.id
    def args = params.vqsr_indel_args ?: ''
    
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" VariantRecalibrator \\
        -R $reference \\
        -V $vcf \\
        --resource:mills,known=false,training=true,truth=true,prior=12.0 $mills \\
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \\
        -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR \\
        -mode INDEL \\
        -O ${prefix}.indel.recal \\
        --tranches-file ${prefix}.indel.tranches \\
        --rscript-file ${prefix}.indel.plots.R \\
        $args
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version | grep GATK | sed 's/GATK v//')
    END_VERSIONS
    """
}

/*
 * PROCESS: Apply VQSR
 */
process APPLY_VQSR {
    tag "$meta.id"
    label 'process_medium'
    container 'broadinstitute/gatk:4.3.0.0'
    
    publishDir "${params.output}/variants/filtered", mode: 'copy'
    
    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta), path(snp_recal), path(snp_tranches)
    tuple val(meta), path(indel_recal), path(indel_tranches)
    path(reference)
    path(reference_dict)
    path(reference_fai)
    
    output:
    tuple val(meta), path("*.vqsr.vcf.gz"), emit: vcf
    tuple val(meta), path("*.vqsr.vcf.gz.tbi"), emit: tbi
    path "versions.yml", emit: versions
    
    when:
    params.apply_vqsr
    
    script:
    def prefix = meta.id
    
    """
    # Apply VQSR for SNPs
    gatk --java-options "-Xmx${task.memory.toGiga()}g" ApplyVQSR \\
        -R $reference \\
        -V $vcf \\
        --recal-file $snp_recal \\
        --tranches-file $snp_tranches \\
        --truth-sensitivity-filter-level 99.0 \\
        --create-output-variant-index true \\
        -mode SNP \\
        -O ${prefix}.snp.vqsr.vcf.gz
    
    # Apply VQSR for INDELs
    gatk --java-options "-Xmx${task.memory.toGiga()}g" ApplyVQSR \\
        -R $reference \\
        -V ${prefix}.snp.vqsr.vcf.gz \\
        --recal-file $indel_recal \\
        --tranches-file $indel_tranches \\
        --truth-sensitivity-filter-level 99.0 \\
        --create-output-variant-index true \\
        -mode INDEL \\
        -O ${prefix}.vqsr.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version | grep GATK | sed 's/GATK v//')
    END_VERSIONS
    """
}

/*
 * PROCESS: Collect variant calling metrics
 */
process VARIANT_METRICS {
    tag "$meta.id"
    label 'process_medium'
    container 'broadinstitute/gatk:4.3.0.0'
    
    publishDir "${params.output}/variants/metrics", mode: 'copy'
    
    input:
    tuple val(meta), path(vcf), path(tbi)
    path(reference)
    path(reference_dict)
    path(reference_fai)
    path(dbsnp)
    path(dbsnp_idx)
    
    output:
    tuple val(meta), path("*.variant_calling_*"), emit: metrics
    tuple val(meta), path("*.variant_calling_metrics.json"), emit: json
    path "versions.yml", emit: versions
    
    script:
    def prefix = meta.id
    
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" CollectVariantCallingMetrics \\
        -I $vcf \\
        -O ${prefix}.variant_calling \\
        --DBSNP $dbsnp \\
        -SD $reference_dict \\
        --THREAD_COUNT $task.cpus
    
    # Create a simple JSON summary (in a real pipeline, this would be more comprehensive)
    echo "{" > ${prefix}.variant_calling_metrics.json
    echo "  \\"sample\\": \\"${meta.id}\\"," >> ${prefix}.variant_calling_metrics.json
    echo "  \\"metrics\\": {" >> ${prefix}.variant_calling_metrics.json
    echo "    \\"snps\\": \\"\$(grep -m 1 'TOTAL_SNPS' ${prefix}.variant_calling_summary_metrics | cut -f 2)\\"," >> ${prefix}.variant_calling_metrics.json
    echo "    \\"indels\\": \\"\$(grep -m 1 'TOTAL_INDELS' ${prefix}.variant_calling_summary_metrics | cut -f 2)\\"," >> ${prefix}.variant_calling_metrics.json
    echo "    \\"pass_filter\\": \\"\$(grep -m 1 'TOTAL_FILTERED' ${prefix}.variant_calling_summary_metrics | cut -f 2)\\"" >> ${prefix}.variant_calling_metrics.json
    echo "  }" >> ${prefix}.variant_calling_metrics.json
    echo "}" >> ${prefix}.variant_calling_metrics.json
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version | grep GATK | sed 's/GATK v//')
    END_VERSIONS
    """
}

// Main workflow for variant calling
workflow VARIANT_CALLING {
    take:
    bam_bai            // [meta, bam, bai]
    reference          // reference fasta
    reference_dict     // reference dictionary
    reference_fai      // reference fasta index
    dbsnp              // dbSNP VCF
    dbsnp_idx          // dbSNP VCF index
    regions            // BED file with regions of interest (optional)
    hapmap             // HapMap resource (for VQSR)
    hapmap_idx         // HapMap index
    omni               // 1000G Omni resource (for VQSR)
    omni_idx           // 1000G Omni index
    thousandG          // 1000G resource (for VQSR)
    thousandG_idx      // 1000G index
    mills              // Mills & 1000G Gold Standard Indels (for VQSR)
    mills_idx          // Mills & 1000G index
    
    main:
    ch_versions = Channel.empty()
    
    // Call variants with appropriate caller
    if (params.variant_caller == 'haplotypecaller') {
        HAPLOTYPECALLER(bam_bai, reference, reference_dict, reference_fai, regions)
        ch_vcf_tbi = HAPLOTYPECALLER.out.vcf.join(HAPLOTYPECALLER.out.tbi, by: [0])
        ch_versions = ch_versions.mix(HAPLOTYPECALLER.out.versions)
    } else if (params.variant_caller == 'mutect2') {
        MUTECT2(bam_bai, reference, reference_dict, reference_fai, regions, [])
        ch_vcf_tbi = MUTECT2.out.vcf.join(MUTECT2.out.tbi, by: [0])
        ch_versions = ch_versions.mix(MUTECT2.out.versions)
    }
    
    // Apply filtering based on strategy
    if (params.apply_vqsr) {
        // VQSR approach - separate recalibration for SNPs and INDELs
        VARIANT_RECALIBRATOR_SNP(
            ch_vcf_tbi,
            reference,
            reference_dict,
            reference_fai,
            hapmap,
            hapmap_idx,
            omni,
            omni_idx,
            thousandG,
            thousandG_idx,
            dbsnp,
            dbsnp_idx
        )
        ch_versions = ch_versions.mix(VARIANT_RECALIBRATOR_SNP.out.versions)
        
        VARIANT_RECALIBRATOR_INDEL(
            ch_vcf_tbi,
            reference,
            reference_dict,
            reference_fai,
            mills,
            mills_idx,
            dbsnp,
            dbsnp_idx
        )
        ch_versions = ch_versions.mix(VARIANT_RECALIBRATOR_INDEL.out.versions)
        
        // Apply recalibration
        ch_snp_recal_tranches = VARIANT_RECALIBRATOR_SNP.out.recal.join(VARIANT_RECALIBRATOR_SNP.out.tranches, by: [0])
        ch_indel_recal_tranches = VARIANT_RECALIBRATOR_INDEL.out.recal.join(VARIANT_RECALIBRATOR_INDEL.out.tranches, by: [0])
        
        APPLY_VQSR(
            ch_vcf_tbi,
            ch_snp_recal_tranches,
            ch_indel_recal_tranches,
            reference,
            reference_dict,
            reference_fai
        )
        ch_filtered_vcf_tbi = APPLY_VQSR.out.vcf.join(APPLY_VQSR.out.tbi, by: [0])
        ch_versions = ch_versions.mix(APPLY_VQSR.out.versions)
    } else {
        // Hard filtering approach
        VARIANT_FILTRATION(
            ch_vcf_tbi,
            reference,
            reference_dict,
            reference_fai
        )
        ch_filtered_vcf_tbi = VARIANT_FILTRATION.out.vcf.join(VARIANT_FILTRATION.out.tbi, by: [0])
        ch_versions = ch_versions.mix(VARIANT_FILTRATION.out.versions)
    }
    
    // Collect variant metrics
    VARIANT_METRICS(
        ch_filtered_vcf_tbi,
        reference,
        reference_dict,
        reference_fai,
        dbsnp,
        dbsnp_idx
    )
    ch_versions = ch_versions.mix(VARIANT_METRICS.out.versions)
    
    emit:
    raw_vcf = ch_vcf_tbi
    filtered_vcf = ch_filtered_vcf_tbi
    metrics = VARIANT_METRICS.out.metrics
    metrics_json = VARIANT_METRICS.out.json
    versions = ch_versions
} 