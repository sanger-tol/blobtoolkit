include { MOSDEPTH } from '../../modules/nf-core/modules/mosdepth/main'

workflow COVERAGE_STATS {
    take:
    cram    // channel: [val(meta), path(cram), path(cai)]
    bed     // path/to/bed
    fasta   // path/to/fasta

    main:
    ch_versions = Channel.empty()

    MOSDEPTH(cram, bed, fasta)
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)

    emit:
    global = MOSDEPTH.out.global_txt
    summary = MOSDEPTH.out.summary_txt
    bed = MOSDEPTH.out.per_base_bed
    base = MOSDEPTH.out.per_base_csi
    versions = ch_versions
}