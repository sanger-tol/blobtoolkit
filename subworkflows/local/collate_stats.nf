//
// Collate genome statistics by various window sizes
//

include { BLOBTOOLKIT_COUNTBUSCOS } from '../../modules/local/blobtoolkit/countbuscos'
include { WINDOWSTATS_INPUT       } from '../../modules/local/windowstats_input'
include { GET_WINDOW_STATS        } from '../../modules/local/get_window_stats'


workflow COLLATE_STATS {
    take: 
    busco_table // channel: [ val(meta), path(full_table) ]
    bed         // channel: [ val(meta), path(bed) ]
    freq        // channel: [ val(meta), path(freq) ]
    mononuc     // channel: [ val(meta), path(mononuc) ]
    cov         // channel: [ val(meta), path(regions.bed.gz) ]

    main:
    ch_versions = Channel.empty()


    // Count BUSCO genes in a region
    busco_table
    | groupTuple()
    | set { ch_busco }

    BLOBTOOLKIT_COUNTBUSCOS ( ch_busco, bed )
    ch_versions = ch_versions.mix ( BLOBTOOLKIT_COUNTBUSCOS.out.versions.first() )


    // Combine outputs from Fasta windows, mosdepth, and count BUSCO genes
    WINDOWSTATS_INPUT ( freq, mononuc, cov, BLOBTOOLKIT_COUNTBUSCOS.out.tsv )
    ch_versions = ch_versions.mix ( WINDOWSTATS_INPUT.out.versions.first() )


    // Genome statistics by different window sizes
    GET_WINDOW_STATS ( WINDOWSTATS_INPUT.out.tsv )
    ch_versions = ch_versions.mix ( GET_WINDOW_STATS.out.versions.first() )


    emit:
    window_tsv = GET_WINDOW_STATS.out.tsv // channel: [ val(meta), path(window_stats_tsvs) ]
    versions   = ch_versions              // channel: [ versions.yml ]
}
