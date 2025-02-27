/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-schema'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowBlobtoolkit.initialise(params, log)

// Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta, params.taxdump, params.busco, params.blastp, params.blastx, params.lineage_tax_ids ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.fasta) { ch_fasta = Channel.value([ [ 'id': params.accession ?: file(params.fasta.replace(".gz", "")).baseName ], file(params.fasta) ]) } else { exit 1, 'Genome fasta file must be specified!' }
if (params.taxon) { ch_taxon = Channel.value(params.taxon) } else { exit 1, 'NCBI Taxon ID not specified!' }
if (params.blastp) { ch_blastp = Channel.fromPath(params.blastp).map { tuple(["type": "blastp"], it) } } else { exit 1, 'Diamond BLASTp database must be specified!' }
if (params.blastx) { ch_blastx = Channel.fromPath(params.blastx).map { tuple(["type": "blastx"], it) } } else { exit 1, 'Diamond BLASTx database must be specified!' }
if (params.blastn) { ch_blastn = Channel.fromPath(params.blastn).map { tuple(["type": "blastn"], it) } } else { exit 1, 'BLASTn database not specified!' }
if (params.taxdump) { ch_taxdump = Channel.fromPath(params.taxdump).map { tuple(["type": "taxdump"], it) } } else { exit 1, 'NCBI Taxonomy database not specified!' }
if (params.fetchngs_samplesheet && !params.align) { exit 1, '--align not specified, even though the input samplesheet is a nf-core/fetchngs one - i.e has fastq files!' }

if (params.lineage_tax_ids) { ch_lineage_tax_ids = Channel.fromPath(params.lineage_tax_ids).first() } else { exit 1, 'Mapping BUSCO lineage <-> taxon_ids not specified' }

// Create channel for optional parameters
if (params.busco_lineages) { ch_busco_lin = Channel.value(params.busco_lineages) } else { ch_busco_lin = Channel.value([]) }
if (params.busco) {
    ch_busco_db = Channel.fromPath(params.busco).first().map { tuple([ "type": "busco"], it ) }
} else {
    ch_busco_db = Channel.value([])
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

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

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow BLOBTOOLKIT {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Prepare genome for downstream processing
    //
    PREPARE_GENOME ( ch_fasta )
    ch_versions = ch_versions.mix ( PREPARE_GENOME.out.versions )

    //
    // SUBWORKFLOW: Check samplesheet and create channels for downstream analysis
    //
    INPUT_CHECK (
        ch_input,
        PREPARE_GENOME.out.genome,
        ch_taxon,
        ch_busco_lin,
        ch_lineage_tax_ids,
        ch_blastn,
        ch_blastx,
        ch_blastp,
        ch_busco_db,
        ch_taxdump,
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
        INPUT_CHECK.out.busco_db,
        INPUT_CHECK.out.blastp,
        INPUT_CHECK.out.taxon_id,
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
    // MODULE: Combine different versions.yml
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // SUBWORKFLOW: Finalise and publish the blobdir
    //
    FINALISE_BLOBDIR (
        BLOBTOOLS.out.blobdir,
        CUSTOM_DUMPSOFTWAREVERSIONS.out.yml,
        VIEW.out.summary
    )
    // Don't update ch_versions because it's already been consumed by now


    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowBlobtoolkit.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowBlobtoolkit.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(BUSCO_DIAMOND.out.multiqc.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
