//
// Take provided samplesheet and genome
// check them and create appropriate channels for downstream analysis
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'
include { GUNZIP            } from '../../modules/nf-core/gunzip/main'

workflow INPUT_CHECK {
    take:
    input     // path(samplesheet)
    fasta     // path(fasta)

    main:
    ch_versions = Channel.empty()

    // Uncompress genome fasta file if required
    if (params.fasta.endsWith('.gz')) {
        genome      = GUNZIP ( fasta.map { file -> [ [ id: file.baseName.replaceFirst(/.fa.*/, "") ], file ] } ).gunzip
        ch_versions = ch_versions.mix(GUNZIP.out.versions)
    } else {
        genome = fasta.map { file -> [ [ id: file.baseName ], file ] }
    }

    SAMPLESHEET_CHECK ( input )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_data_channels(it) }
        .set { aln }
    ch_versions = ch_versions.mix(SAMPLESHEET_CHECK.out.versions)

    emit:
    aln                                       // channel: [ val(meta), [ datafile ] ]
    genome                                    // channel: [ val(meta), fasta ]
    versions = ch_versions                    // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ datafile ] ]
def create_data_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id         = row.sample
    meta.datatype   = row.datatype

    def array = []
    array = [ meta, file(row.datafile) ]
    return array
}
