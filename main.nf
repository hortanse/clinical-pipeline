#!/usr/bin/env nextflow

/*
========================================================================================
    Clinical Bioinformatics Pipeline
========================================================================================
    A robust, clinical-grade genomic analysis pipeline for variant detection and
    interpretation in a healthcare setting.
    
    Homepage: https://github.com/hortanse/clinical-pipeline
========================================================================================
*/

nextflow.enable.dsl = 2

// Print pipeline info
log.info """\
         C L I N I C A L  B I O I N F O R M A T I C S  P I P E L I N E
         =============================================================
         input       : ${params.input}
         output      : ${params.output}
         genome      : ${params.genome}
         """
         .stripIndent()

// Include modules
include { FETCH_SAMPLES }        from './modules/ingestion'
include { FASTQC; TRIM_READS }   from './modules/preprocessing'
include { ALIGN_READS }          from './modules/alignment'
include { CALL_VARIANTS }        from './modules/variant_calling'
include { ANNOTATE_VARIANTS }    from './modules/annotation'
include { FILTER_VARIANTS }      from './modules/interpretation'
include { GENERATE_REPORT }      from './modules/reporting'
include { SEND_TO_EHR }          from './modules/ehr_integration'

// Define workflow
workflow {
    // 1. Fetch sample information from LIMS
    samples_ch = FETCH_SAMPLES(params.input)
    
    // 2. Read QC and preprocessing
    fastqc_ch = FASTQC(samples_ch)
    trimmed_reads_ch = TRIM_READS(samples_ch)
    
    // 3. Alignment to reference genome
    aligned_ch = ALIGN_READS(trimmed_reads_ch, params.genome)
    
    // 4. Call variants
    variants_ch = CALL_VARIANTS(aligned_ch)
    
    // 5. Annotate variants
    annotated_ch = ANNOTATE_VARIANTS(variants_ch)
    
    // 6. Filter and interpret variants
    filtered_ch = FILTER_VARIANTS(annotated_ch)
    
    // 7. Generate clinical report
    report_ch = GENERATE_REPORT(filtered_ch)
    
    // 8. Send report to EHR
    if (params.ehr_integration) {
        SEND_TO_EHR(report_ch)
    }
}

// Error handling
workflow.onError {
    log.error "Pipeline execution stopped with error: ${workflow.errorMessage}"
    // Send notification of failure
    if (params.email) {
        sendMail(
            to: params.email,
            subject: "Clinical Pipeline Execution Failure",
            body: "Pipeline failed with error: ${workflow.errorMessage}"
        )
    }
} 