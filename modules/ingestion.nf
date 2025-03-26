#!/usr/bin/env nextflow

/*
========================================================================================
    INGESTION MODULE
========================================================================================
    This module handles the ingestion of raw data and metadata:
    - Samplesheet parsing
    - FASTQ file handling
    - Input channel creation for downstream processes
*/

// Module parameters with defaults
params.input_format = "fastq"  // Options: fastq, bam
params.samplesheet = null
params.input_dir = null

/*
 * PROCESS: Parse samplesheet and create sample metadata
 */
process PARSE_SAMPLESHEET {
    tag "samplesheet_parsing"
    label 'process_low'
    
    publishDir "${params.output}/pipeline_info", mode: 'copy'
    
    input:
    path samplesheet
    
    output:
    path "*.validated.csv", emit: validated_samplesheet
    path "versions.yml", emit: versions
    
    script:
    """
    #!/usr/bin/env python3
    
    import csv
    import os
    
    # Read samplesheet and validate
    with open("${samplesheet}", 'r') as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames
        
        # Ensure required columns exist
        required_columns = ['sample_id', 'fastq_1', 'fastq_2']
        for col in required_columns:
            if col not in header:
                raise ValueError(f"Required column '{col}' missing from samplesheet")
        
        # Create validated samplesheet
        with open("samplesheet.validated.csv", 'w', newline='') as out_f:
            writer = csv.DictWriter(out_f, fieldnames=header)
            writer.writeheader()
            
            # Write validated rows
            for row in reader:
                # Validate sample_id is not empty
                if not row['sample_id']:
                    continue
                
                # Validate fastq_1 exists
                if not row['fastq_1'] or not os.path.exists(row['fastq_1']):
                    continue
                
                # Write valid row
                writer.writerow(row)
    
    # Create versions file
    with open("versions.yml", "w") as versions_file:
        versions_file.write(\"""${task.process}":
        python: \$(python --version | sed 's/Python //g')
    \""")
    """
}

/*
 * PROCESS: Generate Read Group information for alignment
 */
process GENERATE_READGROUPS {
    tag "$meta.id"
    label 'process_low'
    
    publishDir "${params.output}/ingestion/readgroups", mode: 'copy'
    
    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("*.readgroups.json"), emit: readgroups_json
    path "versions.yml", emit: versions
    
    script:
    def prefix = meta.id
    def single_end = meta.single_end ? "true" : "false"
    
    """
    #!/usr/bin/env python3
    
    import json
    import os
    
    # Create simple read group information
    read_group = {
        "id": "${meta.id}",
        "sample": "${meta.id}",
        "library": "lib_${meta.id}",
        "platform": "ILLUMINA",
        "platform_unit": "unit1",
        "center": "center1"
    }
    
    # Generate RG string for BWA
    rg_string = f"@RG\\\\tID:{read_group['id']}\\\\tSM:{read_group['sample']}\\\\tLB:{read_group['library']}\\\\tPL:{read_group['platform']}"
    read_group["rg_string"] = rg_string
    
    # Write read group information
    with open("${prefix}.readgroups.json", "w") as f:
        json.dump(read_group, f, indent=2)
    
    with open("versions.yml", "w") as versions_file:
        versions_file.write(\"""${task.process}":
        python: \$(python --version | sed 's/Python //g')
    \""")
    """
}

// Main workflow for ingestion
workflow INGESTION {
    take:
    samplesheet
    
    main:
    ch_versions = Channel.empty()
    
    // Parse samplesheet
    PARSE_SAMPLESHEET(samplesheet)
    ch_versions = ch_versions.mix(PARSE_SAMPLESHEET.out.versions)
    
    // Read validated samplesheet
    ch_input = Channel
        .fromPath(PARSE_SAMPLESHEET.out.validated_samplesheet)
        .splitCsv(header: true)
        .map { row -> 
            def meta = [
                id: row.sample_id,
                sample_id: row.sample_id,
                single_end: row.fastq_2 ? false : true
            ]
            
            def reads = row.fastq_2 ? [
                file(row.fastq_1), file(row.fastq_2)
            ] : [
                file(row.fastq_1)
            ]
            
            return [meta, reads]
        }
    
    // Generate read groups
    GENERATE_READGROUPS(ch_input)
    ch_versions = ch_versions.mix(GENERATE_READGROUPS.out.versions)
    
    // Combine reads with readgroups
    ch_reads_with_rg = ch_input.join(GENERATE_READGROUPS.out.readgroups_json)
    
    emit:
    reads = ch_reads_with_rg
    readgroups = GENERATE_READGROUPS.out.readgroups_json
    versions = ch_versions
} 