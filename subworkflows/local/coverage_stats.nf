include { MOSDEPTH      } from '../../modules/nf-core/modules/nf-core/mosdepth/main'
include { SAMTOOLS_VIEW } from '../../modules/nf-core/modules/nf-core/samtools/view/main'


workflow COVERAGE_STATS {
    take: 
    cram    // channel: [val(meta), path(cram), path(cai)]
    //bed     // path/to/bed
    fasta   // path/to/fasta
    qname   // path/to/qname

    main:
    ch_versions = Channel.empty()

    // Convert from CRAM to BAM
    SAMTOOLS_VIEW(cram, fasta, qname)
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions)

    // Calculate Coverage 
    ch_csi = SAMTOOLS_VIEW.out.csi
    ch_bam = SAMTOOLS_VIEW.out.bam.join(ch_csi)
    
    MOSDEPTH(ch_bam, [], [])
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)

    emit:
    bam = SAMTOOLS_VIEW.out.bam
    csi = SAMTOOLS_VIEW.out.csi
    global = MOSDEPTH.out.global_txt
    summary = MOSDEPTH.out.summary_txt
    bed = MOSDEPTH.out.per_base_bed
    base = MOSDEPTH.out.per_base_csi
    versions = ch_versions
}