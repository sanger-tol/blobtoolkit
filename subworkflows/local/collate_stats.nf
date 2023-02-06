include { COUNT_BUSCO_GENES    } from '../../modules/local/count_busco_genes'

workflow COLLATE_STATS {
    take: 
    busco_dir       // channel: [val(meta), path(busco_dir)]
    bed             // channel: [val(meta), path(bed)]

    main:
    ch_versions = Channel.empty()

    ch_tsv_path = GrabBuscoFiles(busco_dir).groupTuple(by: [0])

    // Count Busco Genes
    COUNT_BUSCO_GENES(ch_tsv_path, bed)
    ch_versions = ch_versions.mix(COUNT_BUSCO_GENES.out.versions)

    emit:
    count_genes = COUNT_BUSCO_GENES.out.tsv
    versions = ch_versions
}

process GrabBuscoFiles {
    tag "${meta.id}"
    executor 'local'

    input:
    tuple val(meta), path("in")

    output:
    tuple val(meta), path("in/**/full_table.tsv")

    "true"
}