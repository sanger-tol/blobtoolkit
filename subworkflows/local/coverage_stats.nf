include { MOSDEPTH      } from '../../modules/nf-core/mosdepth/main'
include { SAMTOOLS_VIEW } from '../../modules/nf-core/samtools/view/main'
include { COVERAGE_TSV  } from '../../modules/local/coverage_tsv'
include { GUNZIP        } from '../../modules/nf-core/gunzip/main'

workflow COVERAGE_STATS {
    take: 
    cram    // channel: [val(meta), path(cram), path(cai)]
    fasta   // channel: [val(meta), path(fasta)]
    bed     // channel: [val(meta), path(bed)]
    tsv     // channel: [val(meta), path(tsv)]

    main:
    ch_versions = Channel.empty()

    // Convert from CRAM to BAM
    SAMTOOLS_VIEW(cram, fasta.map{it -> it[1]}, [])
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions)

    // BAM Channel
    ch_csi = SAMTOOLS_VIEW.out.csi
    ch_bam = SAMTOOLS_VIEW.out.bam.join(ch_csi)
    
    // Calculate Coverage 
    MOSDEPTH(ch_bam, bed.map{it -> it[1]}, fasta.map{it -> it[1]})
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)

    // Combine Coverage with TSV from count_buscogenes
    ch_mosdepth = GUNZIP(MOSDEPTH.out.regions_bed).gunzip
    COVERAGE_TSV(ch_mosdepth, tsv)
    ch_versions = ch_versions.mix(COVERAGE_TSV.out.versions)

    emit:
    global = MOSDEPTH.out.global_txt
    summary = MOSDEPTH.out.summary_txt
    bed = MOSDEPTH.out.per_base_bed
    base = MOSDEPTH.out.per_base_csi
    coverage = COVERAGE_TSV.out.cov_tsv
    versions = ch_versions
}