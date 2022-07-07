#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Check mandatory parameters
if (params.input && params.fasta) { inputs = [ file(params.input, checkIfExists: true), file(params.fasta) ] }
else if (params.input && params.project) { inputs = [ params.input, params.project ] }
else { exit 1, 'Input not specified. Please include either a samplesheet or Tree of Life organism and project IDs' }

// 
include { INPUT_CHECK   } from './subworkflows/local/input_check'
include { SAMTOOLS_VIEW } from './modules/local/samtools_view'

workflow {
    Channel.of(inputs).set{ch_input}
    INPUT_CHECK ( ch_input )
    ch_fasta = INPUT_CHECK.out.genome.collect()
    SAMTOOLS_VIEW ( INPUT_CHECK.out.aln, ch_fasta )
}
