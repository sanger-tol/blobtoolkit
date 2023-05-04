/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.fasta, params.taxa_file, params.ncbi_taxdump, params.busco_lineages_path, params.diamondblastp_db ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = Channel.fromPath(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.fasta && params.accession) { ch_fasta = Channel.of([ [ 'id': params.accession ], params.fasta ]).collect() } else { exit 1, 'Genome fasta file and accession must be specified!' }
if (params.taxon) { ch_taxon = Channel.of(params.taxon) } else { exit 1, 'NCBI Taxon ID not specified!' }
if (params.diamondblastp_db) { ch_blastp = Channel.fromPath(params.diamondblastp_db) } else { exit 1, 'Diamond BLASTp database location not specified!' }
if (params.blastp_outext) { ch_outext = Channel.of(params.blastp_outext) } else { exit 1, 'Diamond BLASTp output format not specified!' }
if (params.blastp_cols) { ch_cols = Channel.of(params.blastp_cols) } else { exit 1, 'Diamond BLASTp output columns not specified!' }
if (params.ncbi_taxdump) { ch_taxdump = Channel.fromPath(params.ncbi_taxdump) } else { exit 1, 'NCBI Taxonomy database location not specified!' }

// Create channel for optional parameters
if (params.busco_lineages_path) { ch_busco_db = Channel.fromPath(params.busco_lineages_path) } else { ch_busco_db = Channel.empty() }
if (params.yaml && params.accession) { ch_yaml = Channel.of([ [ 'id': params.accession ], params.yaml ]) } else { ch_yaml = Channel.empty() }


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL/NF-CORE SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { INPUT_CHECK    } from '../subworkflows/local/input_check'
include { COVERAGE_STATS } from '../subworkflows/local/coverage_stats'
include { BUSCO_DIAMOND  } from '../subworkflows/local/busco_diamond_blastp'
include { COLLATE_STATS  } from '../subworkflows/local/collate_stats'
include { BLOBTOOLS      } from '../subworkflows/local/blobtools'
include { VIEW           } from '../subworkflows/local/view'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL/NF-CORE MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Created locally in this pipeline
//

include { BLOBTOOLKIT_CONFIG } from '../modules/local/blobtoolkit/config'

//
// MODULE: Installed directly from nf-core/modules
//

include { GUNZIP                      } from '../modules/nf-core/gunzip/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BLOBTOOLKIT {

    ch_versions = Channel.empty()


    //
    // MODULE: Gunzip fasta file if needed
    //
    if ( params.fasta.endsWith('.gz') ) {
        ch_genome   = GUNZIP ( ch_fasta ).gunzip
        ch_versions = ch_versions.mix ( GUNZIP.out.versions.first() )
    } else {
        ch_genome   = ch_fasta
    }


    //
    // SUBWORKFLOW: Check samplesheet and create channels for downstream analysis
    //
    INPUT_CHECK ( ch_input )
    ch_versions = ch_versions.mix ( INPUT_CHECK.out.versions )


    //
    // SUBWORKFLOW: Calculate genome coverage and statistics 
    //
    COVERAGE_STATS ( INPUT_CHECK.out.aln, ch_genome )
    ch_versions = ch_versions.mix ( COVERAGE_STATS.out.versions )


    //
    // SUBWORKFLOW: Run BUSCO using lineages fetched from GOAT, then run diamond_blastp
    //
    if (params.taxa_file) { 
        ch_taxa = Channel.from(params.taxa_file)
        ch_taxon_taxa = ch_genome.combine(ch_taxon).combine(ch_taxa).map { meta, fasta, taxon, taxa -> [ meta, taxon, taxa ] }
    } else { 
        ch_taxon_taxa = ch_genome.combine(ch_taxon).map { meta, fasta, taxon -> [ meta, taxon, [] ] }
    }

    BUSCO_DIAMOND ( ch_genome, ch_taxon_taxa, ch_busco_db, ch_blastp, ch_outext, ch_cols )
    ch_versions = ch_versions.mix ( BUSCO_DIAMOND.out.versions )


    //
    // SUBWORKFLOW: Collate genome statistics by various window sizes
    //
    COLLATE_STATS ( BUSCO_DIAMOND.out.full_table, COVERAGE_STATS.out.bed, COVERAGE_STATS.out.freq, COVERAGE_STATS.out.mononuc, COVERAGE_STATS.out.cov )
    ch_versions = ch_versions.mix ( COLLATE_STATS.out.versions )


    //
    // SUBWORKFLOW: Create BlobTools dataset
    //
    if ( !params.yaml ) {
        BLOBTOOLKIT_CONFIG ( ch_genome )
        ch_config   = BLOBTOOLKIT_CONFIG.out.yaml
        ch_versions = ch_versions.mix ( BLOBTOOLKIT_CONFIG.out.versions.first() )
    } else {
        ch_config   = ch_yaml
    }

    BLOBTOOLS ( ch_config, COLLATE_STATS.out.window_tsv, BUSCO_DIAMOND.out.first_table, BUSCO_DIAMOND.out.blastp_txt.ifEmpty([[],[]]), ch_taxdump )
    ch_versions = ch_versions.mix ( BLOBTOOLS.out.versions )
    

    //
    // SUBWORKFLOW: Generate static images and summary
    //
    VIEW ( BLOBTOOLS.out.blobdir )
    ch_versions = ch_versions.mix(VIEW.out.versions)
    
    
    //
    // MODULE: Combine different versions.yml
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions
        | unique { it.text }
        | collectFile ( name: 'collated_versions.yml' )
    )
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if ( params.email || params.email_on_fail ) {
        NfcoreTemplate.email ( workflow, params, summary_params, projectDir, log )
    }
    NfcoreTemplate.summary ( workflow, params, log )
    if ( params.hook_url ) {
        NfcoreTemplate.IM_notification ( workflow, params, summary_params, projectDir, log )
    }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
