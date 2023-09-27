//
// Calculate genome coverage and statistics
//

include { SAMTOOLS_VIEW  } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'
include { MOSDEPTH       } from '../../modules/nf-core/mosdepth/main'
include { FASTAWINDOWS   } from '../../modules/nf-core/fastawindows/main'
include { CREATE_BED     } from '../../modules/local/create_bed'


workflow COVERAGE_STATS {
    take: 
    input    // channel: [ val(meta), path(aln) ] 
    fasta    // channel: [ val(meta), path(fasta) ]


    main:
    ch_versions = Channel.empty()


    // Create aligned BAM and index CSI channel
    input
    | branch { meta, aln ->
        bam : aln.toString().endsWith("bam") == true
            return [ meta, aln ]
        cram : aln.toString().endsWith("cram") == true
            return [ meta, aln, [] ]
    }
    | set { ch_aln_idx}

    SAMTOOLS_VIEW ( ch_aln_idx.cram, fasta, [] )
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
    FASTAWINDOWS ( fasta )
    ch_versions = ch_versions.mix ( FASTAWINDOWS.out.versions.first() )


    // Create genome windows file in BED format
    CREATE_BED ( FASTAWINDOWS.out.mononuc )
    ch_versions = ch_versions.mix ( CREATE_BED.out.versions.first() )

    
    // Calculate coverage
    ch_bam_csi
    | combine ( CREATE_BED.out.bed )
    | map { meta, bam, csi, meta2, bed -> [ meta, bam, csi, bed ] }
    | set { ch_bam_csi_bed }
    
    MOSDEPTH ( ch_bam_csi_bed, fasta )
    ch_versions = ch_versions.mix ( MOSDEPTH.out.versions.first() )


    // Combining mosdepth regions_bed in single channel
    MOSDEPTH.out.regions_bed
    | combine ( fasta )
    | map { meta, bed, meta2, fasta -> [ meta2, bed ] }
    | groupTuple ()
    | set { ch_coverage }


    // Mosdepth results for MULTIQC
    MOSDEPTH.out.regions_txt
    | ifEmpty ( MOSDEPTH.out.global_txt )
    | set { multiqc }


    emit:
    freq     = FASTAWINDOWS.out.freq       // channel: [ val(meta), path(freq) ]
    mononuc  = FASTAWINDOWS.out.mononuc    // channel: [ val(meta), path(mononuc) ]
    bed      = CREATE_BED.out.bed          // channel: [ val(meta), path(bed) ]
    cov      = ch_coverage                 // channel: [ val(meta), path(regions.bed.gz) ]
    multiqc                                // channel: [ val(meta), path(dist.txt) ]
    versions = ch_versions                 // channel: [ versions.yml ]
}
