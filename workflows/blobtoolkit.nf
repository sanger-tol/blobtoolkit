/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = Channel.fromPath(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.fasta) { ch_fasta = Channel.fromPath(params.fasta) } else { exit 1, 'Genome fasta file not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow BLOBTOOLKIT {

    ch_versions = Channel.empty()
    ch_ncbi_taxdump = Channel.fromPath(params.ncbi_taxdump)
    blastp_db = Channel.fromPath(params.diamondblastp_db)
    blastp_outext = Channel.of(params.blastp_outext)
    blastp_cols = Channel.of(params.blastp_cols)

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK ( ch_input, ch_fasta )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // SUBWORKFLOW: Convert CRAM to BAM and calculate coverage
    //
    ch_cram = INPUT_CHECK.out.aln.map{ it + [ [] ]}
    ch_fasta = INPUT_CHECK.out.genome
    COVERAGE_STATS(ch_cram, ch_fasta)
    ch_versions = ch_versions.mix(COVERAGE_STATS.out.versions)
    
    //
    // SUBWORKFLOW: Run BUSCO using lineages fetched from GOAT, then run diamond_blastp
    //
    BUSCO_DIAMOND (
    ch_fasta,
    blastp_db,
    blastp_outext,
    blastp_cols
    )
    ch_versions = ch_versions.mix(BUSCO_DIAMOND.out.versions)

    //
    // SUBWORKFLOW: Count BUSCO genes   
    //
    COLLATE_STATS(BUSCO_DIAMOND.out.busco_dir, COVERAGE_STATS.out.fw_bed, COVERAGE_STATS.out.regions_bed)
    ch_versions = ch_versions.mix(COLLATE_STATS.out.versions)

    //
    // SUBWORKFLOW: BLOBTOOLS
    //
    BLOBTOOLS(
    ch_fasta,
    COLLATE_STATS.out.window_tsv,
    BUSCO_DIAMOND.out.first_table,
    BUSCO_DIAMOND.out.blastp_txt,
    INPUT_CHECK.out.genome.map{ meta, fasta -> fasta.baseName },
    ch_ncbi_taxdump
    )
    ch_versions = ch_versions.mix(BLOBTOOLS.out.versions)

    //
    // MODULE: Combine different versions.yml
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
    ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
