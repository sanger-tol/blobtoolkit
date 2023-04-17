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
    // Channel: [meta, cram, cai, meta2, fasta]
    input_sam = cram.combine(fasta)
    SAMTOOLS_VIEW( 
        input_sam.map{ meta, cram, cai, meta2, fasta -> [ meta, cram, cai ] },
        input_sam.map{ meta, cram, cai, meta2, fasta -> fasta },
        []
    )
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
    // Channel: [meta, bam, csi, meta2, bed, meta3, fasta]
    bam_bed = ch_bam.combine(ch_bed)
    input_depth = bam_bed.combine(fasta)
    MOSDEPTH(
        input_depth.map{ meta, bam, csi, meta2, bed, meta3, fasta -> [ meta, bam, csi, bed ] },
        input_depth.map{ meta, bam, csi, meta2, bed, meta3, fasta -> [ meta3, fasta ] }
    )
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)

    emit:
    global = MOSDEPTH.out.global_txt
    summary = MOSDEPTH.out.summary_txt
    regions_bed = MOSDEPTH.out.regions_bed
    regions_csi = MOSDEPTH.out.regions_csi
    fw_bed = ch_bed
    versions = ch_versions
}