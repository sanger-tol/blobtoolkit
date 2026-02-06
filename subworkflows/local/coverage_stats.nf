//
// Calculate genome coverage and statistics
//

include { SAMTOOLS_VIEW  } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'
include { BLOBTK_DEPTH   } from '../../modules/local/blobtk/depth'
include { FASTAWINDOWS   } from '../../modules/nf-core/fastawindows/main'
include { PIGZ_COMPRESS  } from '../../modules/nf-core/pigz/compress/main'
include { CREATE_BED     } from '../../modules/local/create_bed'


workflow COVERAGE_STATS {
    take:
    ch_reads    // channel: [ val(meta), path(aln) ]
    ch_fasta    // channel: [ val(meta), path(fasta) ]


    main:
    ch_versions = channel.empty()


    // Create aligned BAM and index CSI channel
    ch_reads
    | branch { meta, aln ->
        bam : aln.name.endsWith("bam")
            return [ meta, aln ]
        cram : aln.name.endsWith("cram")
            return [ meta, aln, [] ]
    }
    | set { ch_aln_idx}

    SAMTOOLS_VIEW ( ch_aln_idx.cram, ch_fasta, [] )
    ch_versions = ch_versions.mix ( SAMTOOLS_VIEW.out.versions.first() )

    SAMTOOLS_VIEW.out.bam
    | join ( SAMTOOLS_VIEW.out.csi )
    | set { ch_view }

    SAMTOOLS_INDEX ( ch_aln_idx.bam )
    ch_versions = ch_versions.mix ( SAMTOOLS_INDEX.out.versions.first() )

    ch_aln_idx.bam
    | join ( SAMTOOLS_INDEX.out.csi )
    | set { ch_index }

    ch_view
    | mix ( ch_index )
    | set { ch_bam_csi }


    // Calculate genome statistics
    FASTAWINDOWS ( ch_fasta )
    ch_versions = ch_versions.mix ( FASTAWINDOWS.out.versions.first() )


    // Compress the TSV files
    PIGZ_COMPRESS (
        FASTAWINDOWS.out.mononuc
        | mix ( FASTAWINDOWS.out.dinuc )
        | mix ( FASTAWINDOWS.out.trinuc )
        | mix ( FASTAWINDOWS.out.tetranuc )
        | mix ( FASTAWINDOWS.out.freq )
    )
    ch_versions = ch_versions.mix ( PIGZ_COMPRESS.out.versions.first() )


    // Create genome windows file in BED format
    CREATE_BED ( FASTAWINDOWS.out.mononuc )
    ch_versions = ch_versions.mix ( CREATE_BED.out.versions.first() )


    // Calculate coverage
    BLOBTK_DEPTH ( ch_bam_csi )
    ch_versions = ch_versions.mix ( BLOBTK_DEPTH.out.versions.first() )


    // Combining  regions_bed in single channel
    BLOBTK_DEPTH.out.bed
    | combine ( ch_fasta )
    | map { _meta, bed, meta2, _fasta -> [ meta2, bed ] }
    | groupTuple ()
    | set { ch_coverage }


    emit:
    freq     = FASTAWINDOWS.out.freq       // channel: [ val(meta), path(freq) ]
    mononuc  = FASTAWINDOWS.out.mononuc    // channel: [ val(meta), path(mononuc) ]
    bed      = CREATE_BED.out.bed          // channel: [ val(meta), path(bed) ]
    cov      = ch_coverage                 // channel: [ val(meta), path(regions.bed.gz) ]
    versions = ch_versions                 // channel: [ versions.yml ]
}
