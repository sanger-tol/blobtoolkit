//
// Calculate genome coverage and statistics
//

include { SAMTOOLS_VIEW  } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'
include { BLOBTK_DEPTH   } from '../../modules/nf-core/blobtk/depth'
include { FASTAWINDOWS   } from '../../modules/nf-core/fastawindows/main'
include { PIGZ_COMPRESS  } from '../../modules/nf-core/pigz/compress/main'
include { CREATE_BED     } from '../../modules/local/create_bed'


workflow COVERAGE_STATS {
    take:
    ch_reads    // channel: [ val(meta), path(aln) ]
    ch_fasta    // channel: [ val(meta), path(fasta) ]


    main:
    ch_versions = channel.empty()


    //
    // MODULE: CREATE ALIGNED BAM AND INDEX CSI CHANNEL
    //
    ch_aln_idx = ch_reads
        .branch { meta, aln ->
            bam : aln.name.endsWith("bam")
                return [ meta, aln ]
            cram : aln.name.endsWith("cram")
                return [ meta, aln, [] ]
        }
    ch_fa_fai = ch_fasta
        .map { meta, fa -> [ meta, fa, [] ] }

    SAMTOOLS_VIEW ( ch_aln_idx.cram, ch_fa_fai, [], false )
    ch_view = SAMTOOLS_VIEW.out.bam.join(SAMTOOLS_VIEW.out.csi)


    //
    // MODULE: INDEX ALIGNED BAM
    //
    SAMTOOLS_INDEX ( ch_aln_idx.bam )
    ch_index = ch_aln_idx.bam.join(SAMTOOLS_INDEX.out.index)
    ch_bam_csi = ch_view.mix(ch_index)


    //
    // MODULE: CALCULATE GENOME STATISTICS
    //
    FASTAWINDOWS ( ch_fasta )


    //
    // MODULE: COMPRESS TSV FILES
    //
    PIGZ_COMPRESS (
        FASTAWINDOWS.out.mononuc.mix(FASTAWINDOWS.out.dinuc).mix(FASTAWINDOWS.out.trinuc).mix(FASTAWINDOWS.out.tetranuc).mix(FASTAWINDOWS.out.freq)
    )


    //
    // MODULE: CREATE GENOME WINDOWS FILE IN BED FORMAT
    //
    CREATE_BED ( FASTAWINDOWS.out.mononuc )
    ch_versions = ch_versions.mix ( CREATE_BED.out.versions.first() )


    //
    // MODULE: CALCULATE COVERAGE
    //
    BLOBTK_DEPTH ( ch_bam_csi )


    //
    // LOGIC: COMBINE THE BLOBTK_DEPTH BED FILES INTO 1 CHANNEL PER meta
    //
    ch_coverage = BLOBTK_DEPTH.out.bed
        .combine(ch_fasta)
        .map { _meta, bed, meta2, _fasta -> [ meta2, bed ] }
        .groupTuple()


    emit:
    freq     = FASTAWINDOWS.out.freq       // channel: [ val(meta), path(freq) ]
    mononuc  = FASTAWINDOWS.out.mononuc    // channel: [ val(meta), path(mononuc) ]
    bed      = CREATE_BED.out.bed          // channel: [ val(meta), path(bed) ]
    cov      = ch_coverage                 // channel: [ val(meta), path(regions.bed.gz) ]
    versions = ch_versions                 // channel: [ versions.yml ]
}
