/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { PREPARE_GENOME     } from '../subworkflows/local/prepare_genome'
include { MINIMAP2_ALIGNMENT } from '../subworkflows/local/minimap_alignment'
include { INPUT_CHECK        } from '../subworkflows/local/input_check'
include { COVERAGE_STATS     } from '../subworkflows/local/coverage_stats'
include { BUSCO_DIAMOND      } from '../subworkflows/local/busco_diamond_blastp'
include { RUN_BLASTX         } from '../subworkflows/local/run_blastx'
include { RUN_BLASTN         } from '../subworkflows/local/run_blastn'
include { COLLATE_STATS      } from '../subworkflows/local/collate_stats'
include { BLOBTOOLS          } from '../subworkflows/local/blobtools'
include { VIEW               } from '../subworkflows/local/view'
include { FINALISE_BLOBDIR   } from '../subworkflows/local/finalise_blobdir'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_blobtoolkit_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BLOBTOOLKIT {

    take:
    ch_fasta
    ch_databases
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    //
    // SUBWORKFLOW: Prepare genome for downstream processing
    //
    PREPARE_GENOME ( ch_fasta )
    ch_versions = ch_versions.mix ( PREPARE_GENOME.out.versions )

    //
    // SUBWORKFLOW: Check samplesheet and create channels for downstream analysis
    //
    INPUT_CHECK (
        params.input,
        PREPARE_GENOME.out.genome,
        params.taxon,
        Channel.value(params.busco_lineages ?: []),
        params.lineage_tax_ids,
        ch_databases,
    )
    ch_versions = ch_versions.mix ( INPUT_CHECK.out.versions )

    //
    // SUBWORKFLOW: Optional read alignment
    //
    if ( params.align ) {
        MINIMAP2_ALIGNMENT ( INPUT_CHECK.out.reads, PREPARE_GENOME.out.genome )
        ch_versions = ch_versions.mix ( MINIMAP2_ALIGNMENT.out.versions )
        ch_aligned = MINIMAP2_ALIGNMENT.out.aln
    } else {
        ch_aligned = INPUT_CHECK.out.reads
    }

    //
    // SUBWORKFLOW: Calculate genome coverage and statistics
    //
    COVERAGE_STATS ( ch_aligned, PREPARE_GENOME.out.genome )
    ch_versions = ch_versions.mix ( COVERAGE_STATS.out.versions )

    //
    // SUBWORKFLOW: Run BUSCO using lineages fetched from GoaT, then run diamond_blastp
    //
    BUSCO_DIAMOND (
        PREPARE_GENOME.out.genome,
        INPUT_CHECK.out.busco_lineages,
        INPUT_CHECK.out.busco_db.first(),
        INPUT_CHECK.out.blastp.first(),
        INPUT_CHECK.out.taxon_id,
        INPUT_CHECK.out.busco_output,
    )
    ch_versions = ch_versions.mix ( BUSCO_DIAMOND.out.versions )

    //
    // SUBWORKFLOW: Diamond blastx search of assembly contigs against the UniProt reference proteomes
    //
    RUN_BLASTX (
        PREPARE_GENOME.out.genome,
        BUSCO_DIAMOND.out.first_table,
        INPUT_CHECK.out.blastx,
        INPUT_CHECK.out.taxon_id,
    )
    ch_versions = ch_versions.mix ( RUN_BLASTX.out.versions )


    //
    // SUBWORKFLOW: Run blastn search on sequences that had no blastx hits
    //
    RUN_BLASTN (
        RUN_BLASTX.out.blastx_out,
        PREPARE_GENOME.out.genome,
        INPUT_CHECK.out.blastn,
        INPUT_CHECK.out.taxon_id,
    )

    //
    // SUBWORKFLOW: Collate genome statistics by various window sizes
    //
    COLLATE_STATS (
        BUSCO_DIAMOND.out.all_tables,
        COVERAGE_STATS.out.bed,
        COVERAGE_STATS.out.freq,
        COVERAGE_STATS.out.mononuc,
        COVERAGE_STATS.out.cov
    )
    ch_versions = ch_versions.mix ( COLLATE_STATS.out.versions )

    //
    // SUBWORKFLOW: Create BlobTools dataset
    //
    BLOBTOOLS (
        INPUT_CHECK.out.config,
        INPUT_CHECK.out.synonyms_tsv.ifEmpty([[],[]]),
        INPUT_CHECK.out.categories_tsv.ifEmpty([[],[]]),
        COLLATE_STATS.out.window_tsv,
        BUSCO_DIAMOND.out.all_tables,
        BUSCO_DIAMOND.out.blastp_txt.ifEmpty([[],[]]),
        RUN_BLASTX.out.blastx_out.ifEmpty([[],[]]),
        RUN_BLASTN.out.blastn_out.ifEmpty([[],[]]),
        INPUT_CHECK.out.taxdump
    )
    ch_versions = ch_versions.mix ( BLOBTOOLS.out.versions )

    //
    // SUBWORKFLOW: Generate summary and static images
    //
    VIEW ( BLOBTOOLS.out.blobdir )
    ch_versions = ch_versions.mix(VIEW.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info/blobtoolkit",
            name:  'blobtoolkit_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // SUBWORKFLOW: Finalise and publish the blobdir
    //
    FINALISE_BLOBDIR (
        BLOBTOOLS.out.blobdir,
        ch_collated_versions,
        VIEW.out.summary
    )
    // Don't update ch_versions because it's already been consumed by now


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    ch_multiqc_files                      = ch_multiqc_files.mix(BUSCO_DIAMOND.out.multiqc.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
