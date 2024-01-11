//
// Check input samplesheet and get aligned read channels
//

include { SAMPLESHEET_CHECK         } from '../../modules/local/samplesheet_check'
include { FETCHNGSSAMPLESHEET_CHECK } from '../../modules/local/fetchngssamplesheet_check'
include { BLOBTOOLKIT_CONFIG        } from '../../modules/local/blobtoolkit/config'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    fasta       // channel: [ meta, path(fasta) ]
    yaml        // channel: [ meta, path(config ]

    main:
    ch_versions = Channel.empty()

    if ( params.fetchngs_samplesheet ) {
        FETCHNGSSAMPLESHEET_CHECK ( samplesheet )
            .csv
            .splitCsv ( header:true, sep:',' )
            .map { create_data_channels_from_fetchngs(it) }
            .set { aln }
        ch_versions = ch_versions.mix ( FETCHNGSSAMPLESHEET_CHECK.out.versions.first() )

    } else {
        SAMPLESHEET_CHECK ( samplesheet )
            .csv
            .splitCsv ( header:true, sep:',' )
            .map { create_data_channels(it) }
            .set { aln }
        ch_versions = ch_versions.mix ( SAMPLESHEET_CHECK.out.versions.first() )
    }


    if ( !params.yaml ) {
        aln
        | map { meta, data -> meta.id.split("_")[0..-2].join("_") }
        | combine ( fasta )
        | map { sample, meta, fasta -> [ meta, sample ] }
        | groupTuple()
        | set { reads }

        BLOBTOOLKIT_CONFIG ( reads, fasta )
        ch_versions = ch_versions.mix ( BLOBTOOLKIT_CONFIG.out.versions.first() )
        ch_config = BLOBTOOLKIT_CONFIG.out.yaml
    } else {
        ch_config   = yaml
    }

    emit:
    aln                                     // channel: [ val(meta), path(datafile) ]
    config = ch_config                      // channel: [ val(meta), path(yaml) ]
    versions = ch_versions                  // channel: [ versions.yml ]
}

// Function to get list of [ meta, datafile ]
def create_data_channels(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.datatype   = row.datatype


    // add path(s) of the read file(s) to the meta map
    def data_meta = []

    if ( !params.align && (row.datafile.endsWith(".fastq") || row.datafile.endsWith(".fastq.gz")) ) {
        exit 1, "ERROR: Please check input samplesheet and pipeline parameters -> Data file is in FastQ format but --align is not set!\n${row.datafile}"
    }

    if ( !file(row.datafile).exists() ) {
        exit 1, "ERROR: Please check input samplesheet -> Data file does not exist!\n${row.datafile}"
    } else {
        data_meta = [ meta, file(row.datafile) ]
    }

    return data_meta
}

// Function to get list of [ meta, datafile ]
def create_data_channels_from_fetchngs(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.run_accession

    switch (row.instrument_platform) {
        case "ILLUMINA":
            meta.datatype = (row.library_strategy == "Hi-C" ? "hic" : "illumina")
            break
        case "OXFORD_NANOPORE":
            meta.datatype = "ont"
            break
        case "PACBIO_SMRT":
            meta.datatype = (row.instrument_model == "Sequel" ? "pacbio_clr" : "pacbio")
            break
        default:
            meta.datatype = "illumina"
    }


    // add path(s) of the read file(s) to the meta map
    def data_meta = []

    if ( !file(row.fastq_1).exists() ) {
        exit 1, "ERROR: Please check input samplesheet -> Data file does not exist!\n${row.fastq_1}"
    } else {
        data_meta = [ meta, file(row.fastq_1) ]
    }

    return data_meta
}

