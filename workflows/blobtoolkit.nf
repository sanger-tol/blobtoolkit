/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowBlobtoolkit.initialise(params, log)

// Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta, params.taxa_file, params.taxdump, params.busco, params.blastp, params.blastx ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.fasta && params.accession) { ch_fasta = Channel.of([ [ 'id': params.accession ], params.fasta ]).first() } else { exit 1, 'Genome fasta file and accession must be specified!' }
if (params.taxon) { ch_taxon = Channel.of(params.taxon) } else { exit 1, 'NCBI Taxon ID not specified!' }
if (params.blastp) { ch_blastp = file(params.blastp) } else { exit 1, 'Diamond BLASTp database not specified!' }
if (params.blastx) { ch_blastx = file(params.blastx) } else { exit 1, 'Diamond BLASTx database not specified!' }
if (params.taxdump) { ch_taxdump = file(params.taxdump) } else { exit 1, 'NCBI Taxonomy database not specified!' }

// Create channel for optional parameters
if (params.busco) { ch_busco_db = Channel.fromPath(params.busco) } else { ch_busco_db = Channel.empty() }
if (params.yaml && params.accession) { ch_yaml = Channel.of([ [ 'id': params.accession ], params.yaml ]) } else { ch_yaml = Channel.empty() }

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
// MODULE: Loaded from modules/local/
//
include { BLOBTOOLKIT_CONFIG } from '../modules/local/blobtoolkit/config'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { PREPARE_GENOME     } from '../subworkflows/local/prepare_genome'
include { MINIMAP2_ALIGNMENT } from '../subworkflows/local/minimap_alignment'
include { INPUT_CHECK        } from '../subworkflows/local/input_check'
include { COVERAGE_STATS     } from '../subworkflows/local/coverage_stats'
include { BUSCO_DIAMOND      } from '../subworkflows/local/busco_diamond_blastp'
include { RUN_BLASTX         } from '../subworkflows/local/run_blastx'
include { COLLATE_STATS      } from '../subworkflows/local/collate_stats'
include { BLOBTOOLS          } from '../subworkflows/local/blobtools'
include { VIEW               } from '../subworkflows/local/view'

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
    INPUT_CHECK ( ch_input )
    ch_versions = ch_versions.mix ( INPUT_CHECK.out.versions )

    // 
    // SUBWORKFLOW: Optional read alignment
    //
    if ( params.align ) {
        MINIMAP2_ALIGNMENT ( INPUT_CHECK.out.aln, PREPARE_GENOME.out.genome )
        ch_versions = ch_versions.mix ( MINIMAP2_ALIGNMENT.out.versions )
        ch_aligned = MINIMAP2_ALIGNMENT.out.aln
    } else {
        ch_aligned = INPUT_CHECK.out.aln
    }

    //
    // SUBWORKFLOW: Calculate genome coverage and statistics 
    //
    COVERAGE_STATS ( ch_aligned, PREPARE_GENOME.out.genome )
    ch_versions = ch_versions.mix ( COVERAGE_STATS.out.versions )

    //
    // SUBWORKFLOW: Run BUSCO using lineages fetched from GOAT, then run diamond_blastp
    //
    if (params.taxa_file) { 
        ch_taxa = Channel.from(params.taxa_file)
        ch_taxon_taxa = ch_fasta.combine(ch_taxon).combine(ch_taxa).map { meta, fasta, taxon, taxa -> [ meta, taxon, taxa ] }
    } else { 
        ch_taxon_taxa = ch_fasta.combine(ch_taxon).map { meta, fasta, taxon -> [ meta, taxon, [] ] }
    }

    BUSCO_DIAMOND ( 
        PREPARE_GENOME.out.genome, 
        ch_taxon_taxa, 
        ch_busco_db, 
        ch_blastp, 
        params.blastp_outext, 
        params.blastp_cols 
    )
    ch_versions = ch_versions.mix ( BUSCO_DIAMOND.out.versions )
    
    //
    // SUBWORKFLOW: Diamond blastx search of assembly contigs against the UniProt reference proteomes
    RUN_BLASTX ( 
        ch_genome,
        BUSCO_DIAMOND.out.first_table,
        ch_blastx,
        params.blastx_outext,
        params.blastx_cols
    )
    ch_versions = ch_versions.mix ( RUN_BLASTX.out.versions )

    //
    // SUBWORKFLOW: Collate genome statistics by various window sizes
    //
    COLLATE_STATS ( 
        BUSCO_DIAMOND.out.full_table, 
        COVERAGE_STATS.out.bed, 
        COVERAGE_STATS.out.freq, 
        COVERAGE_STATS.out.mononuc, 
        COVERAGE_STATS.out.cov 
    )
    ch_versions = ch_versions.mix ( COLLATE_STATS.out.versions )

    //
    // SUBWORKFLOW: Create BlobTools dataset
    //
    if ( !params.yaml ) {
        BLOBTOOLKIT_CONFIG ( PREPARE_GENOME.out.genome )
        ch_config   = BLOBTOOLKIT_CONFIG.out.yaml
        ch_versions = ch_versions.mix ( BLOBTOOLKIT_CONFIG.out.versions.first() )
    } else {
        ch_config   = ch_yaml
    }

    BLOBTOOLS ( 
        ch_config, 
        COLLATE_STATS.out.window_tsv, 
        BUSCO_DIAMOND.out.first_table, 
        BUSCO_DIAMOND.out.blastp_txt.ifEmpty([[],[]]), 
        RUN_BLASTX.out.blastx_out.ifEmpty([[],[]]),
        ch_taxdump
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
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowBlobtoolkit.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowBlobtoolkit.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(BUSCO_DIAMOND.out.multiqc.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(COVERAGE_STATS.out.multiqc.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
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
