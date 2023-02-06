include { MOSDEPTH      } from '../../modules/nf-core/mosdepth/main'
include { SAMTOOLS_VIEW } from '../../modules/nf-core/samtools/view/main'
include { FASTAWINDOWS  } from '../../modules/nf-core/fastawindows/main'
include { CREATE_BED    } from '../../modules/local/create_bed'

workflow COVERAGE_STATS {
    take: 
    cram    // channel: [val(meta), path(cram), path(cai)]
    fasta   // channel: [val(meta), path(fasta)]

    main:
    ch_versions = Channel.empty()

    // Convert from CRAM to BAM
    SAMTOOLS_VIEW(cram, fasta.map{it -> it[1]}, [])
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions)

    // Generate BED File
    FASTAWINDOWS(fasta)
    ch_versions = ch_versions.mix(FASTAWINDOWS.out.versions)

    CREATE_BED(FASTAWINDOWS.out.mononuc)
    ch_versions = ch_versions.mix(CREATE_BED.out.versions)
    
    ch_bed = CREATE_BED.out.bed

    // BAM Channel
    ch_csi = SAMTOOLS_VIEW.out.csi
    ch_bam = SAMTOOLS_VIEW.out.bam.join(ch_csi)
    
    // Calculate Coverage (need to remove `meta` from the `ch_bed` and `fasta` channels)
    MOSDEPTH(ch_bam, ch_bed.map{it -> it[1]}, fasta.map{it -> it[1]})
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)

    emit:
    global = MOSDEPTH.out.global_txt
    summary = MOSDEPTH.out.summary_txt
    base_bed = MOSDEPTH.out.per_base_bed
    base_csi = MOSDEPTH.out.per_base_csi
    bed = ch_bed
    versions = ch_versions
}