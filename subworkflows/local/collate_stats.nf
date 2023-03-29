include { COUNT_BUSCO_GENES    } from '../../modules/local/count_busco_genes'
include { GUNZIP               } from '../../modules/nf-core/gunzip/main'
include { COVERAGE_TSV         } from '../../modules/local/coverage_tsv'
include { GET_WINDOW_STATS     } from '../../modules/local/get_window_stats'


workflow COLLATE_STATS {
    take: 
    busco_dir       // channel: [val(meta), path(busco_dir)]
    bed             // channel: [val(meta), path(bed)]
    regions_bed     // channel: [val(meta), path(regions_bed)]

    main:
    ch_versions = Channel.empty()

    // Extract Busco's full_table.tsv files from the directories
    ch_tsv_path = busco_dir.map {
        meta, busco_dir -> [meta, file("${busco_dir}/**/full_table.tsv")[0]]
    }.groupTuple(by: [0])

    // Count Busco Genes
    COUNT_BUSCO_GENES(ch_tsv_path, bed)
    ch_versions = ch_versions.mix(COUNT_BUSCO_GENES.out.versions)

    // Combine output TSV from mosdepth and count_busco_genes
    COVERAGE_TSV(GUNZIP(regions_bed).gunzip, COUNT_BUSCO_GENES.out.tsv)
    ch_versions = ch_versions.mix(COVERAGE_TSV.out.versions)

    GET_WINDOW_STATS(COVERAGE_TSV.out.cov_tsv)
    ch_versions = ch_versions.mix(GET_WINDOW_STATS.out.versions)

    emit:
    count_genes = COUNT_BUSCO_GENES.out.tsv
    cov_tsv = COVERAGE_TSV.out.cov_tsv
    window_tsv = GET_WINDOW_STATS.out.tsv
    versions = ch_versions
}
