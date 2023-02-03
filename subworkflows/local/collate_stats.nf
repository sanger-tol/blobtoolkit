include { COUNT_BUSCO_GENES    } from '../../modules/local/count_busco_genes'

workflow COLLATE_STATS {
    take: 
    tsv    // channel: [val(meta), path(tsv)]
    bed    // channel: [val(meta), path(bed)]

    main:
    ch_versions = Channel.empty()

    ch_tsv_path = GrabFiles(tsv).groupTuple(by: [0])

    // Count Busco Genes
    COUNT_BUSCO_GENES(ch_tsv_path, bed)
    ch_versions = ch_versions.mix(COUNT_BUSCO_GENES.out.versions)

    emit:
    count_genes = COUNT_BUSCO_GENES.out.tsv
    versions = ch_versions
}

process GrabFiles {
    tag "${meta.id}"
    executor 'local'

    input:
    tuple val(meta), path("in")

    output:
    tuple val(meta), path("in/**/full_table.tsv")

    "true"
}