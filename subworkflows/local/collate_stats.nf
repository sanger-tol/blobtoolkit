include { COUNT_BUSCO_GENES    } from '../../modules/local/count_busco_genes'

workflow COLLATE_STATS {
    take: 
    busco_dir       // channel: [val(meta), path(busco_dir)]
    bed             // channel: [val(meta), path(bed)]

    main:
    ch_versions = Channel.empty()

    // Extract Busco's full_table.tsv files from the directories
    ch_tsv_path = busco_dir.map {
        meta, busco_dir -> [meta, file("${busco_dir}/**/full_table.tsv")[0]]
    }.groupTuple(by: [0])

    // Count Busco Genes
    COUNT_BUSCO_GENES(ch_tsv_path, bed)
    ch_versions = ch_versions.mix(COUNT_BUSCO_GENES.out.versions)

    emit:
    count_genes = COUNT_BUSCO_GENES.out.tsv
    versions = ch_versions
}
