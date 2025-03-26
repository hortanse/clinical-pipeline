/*
========================================================================================
    Variant Calling Module
========================================================================================
    This module handles variant calling using GATK HaplotypeCaller and variant filtering
*/

process CALL_VARIANTS {
    tag "${meta.id}"
    label 'process_high'
    
    input:
    tuple val(meta), path(bam), path(bai)
    path genome
    path regions
    
    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    path "versions.yml", emit: versions
    
    script:
    def prefix = "${meta.id}"
    def avail_mem = task.memory ? "-Xmx${task.memory.toGiga()}g" : ""
    def intervals_command = regions ? "--intervals ${regions}" : ""
    
    """
    # Call variants using HaplotypeCaller
    gatk --java-options "${avail_mem}" HaplotypeCaller \\
        -R ${genome}/genome.fa \\
        -I ${bam} \\
        ${intervals_command} \\
        -O ${prefix}.vcf.gz \\
        --native-pair-hmm-threads ${task.cpus} \\
        --create-output-variant-index true
    
    # Record software versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version | grep -e "^GATK" | sed 's/GATK //g')
    END_VERSIONS
    """
}

process VARIANT_FILTRATION {
    tag "${meta.id}"
    label 'process_medium'
    
    input:
    tuple val(meta), path(vcf), path(tbi)
    path genome
    
    output:
    tuple val(meta), path("*.filtered.vcf.gz"), emit: vcf
    tuple val(meta), path("*.filtered.vcf.gz.tbi"), emit: tbi
    path "versions.yml", emit: versions
    
    script:
    def prefix = "${meta.id}"
    def avail_mem = task.memory ? "-Xmx${task.memory.toGiga()}g" : ""
    def intervals_command = regions ? "--intervals ${regions}" : ""
    def tbis = "${prefix}.vcf.gz.tbi"
    def genome = "${genome}/genome.fa"
    
    """
    # Filter variants
    gatk --java-options "${avail_mem}" VariantFiltration \\
        -R ${genome} \\
        -V ${vcf} \\
        -O ${prefix}.filtered.vcf.gz \\
        --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \\
        --filter-name "hard_filter"
    
    # Index filtered VCF
    gatk --java-options "${avail_mem}" IndexFeatureFile \\
        -I ${prefix}.filtered.vcf.gz \\
        -O ${tbis}
    
    # Record software versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version | grep -e "^GATK" | sed 's/GATK //g')
    END_VERSIONS
    """
}

process COMBINE_GVCFs {
    tag "${meta.id}"
    label 'process_high'
    
    input:
    tuple val(meta), path(gvcfs), path(tbis)
    path genome
    
    output:
    tuple val(meta), path("*.combined.g.vcf.gz"), emit: gvcf
    tuple val(meta), path("*.combined.g.vcf.gz.tbi"), emit: tbi
    path "versions.yml", emit: versions
    
    script:
    def prefix = "${meta.id ?: 'combined'}"
    def avail_mem = task.memory ? "-Xmx${task.memory.toGiga()}g" : ""
    def gvcfs_command = gvcfs.collect { "-V $it" }.join(' ')
    
    """
    # Combine GVCFs
    gatk --java-options "${avail_mem}" CombineGVCFs \\
        -R ${genome}/genome.fa \\
        ${gvcfs_command} \\
        -O ${prefix}.combined.g.vcf.gz
    
    # Record software versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version | grep -e "^GATK" | sed 's/GATK //g')
    END_VERSIONS
    """
}

process GENOTYPE_GVCF {
    tag "${meta.id}"
    label 'process_high'
    
    input:
    tuple val(meta), path(gvcf), path(tbi)
    path genome
    
    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    path "versions.yml", emit: versions
    
    script:
    def prefix = "${meta.id}"
    def avail_mem = task.memory ? "-Xmx${task.memory.toGiga()}g" : ""
    
    """
    # Joint genotyping
    gatk --java-options "${avail_mem}" GenotypeGVCFs \\
        -R ${genome}/genome.fa \\
        -V ${gvcf} \\
        -O ${prefix}.vcf.gz
    
    # Record software versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version | grep -e "^GATK" | sed 's/GATK //g')
    END_VERSIONS
    """
} 