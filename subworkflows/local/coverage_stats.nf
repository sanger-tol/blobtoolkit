//
// Calculate genome coverage and statistics
//

include { SAMTOOLS_VIEW } from '../../modules/nf-core/samtools/view/main'
include { MOSDEPTH      } from '../../modules/nf-core/mosdepth/main'
include { FASTAWINDOWS  } from '../../modules/nf-core/fastawindows/main'
include { CREATE_BED    } from '../../modules/local/create_bed'


workflow COVERAGE_STATS {
    take: 
    cram    // channel: [ val(meta), path(cram) ] 
    fasta   // channel: [ val(meta), path(fasta) ]


    main:
    ch_versions = Channel.empty()


    // Convert from CRAM to BAM
    cram
    | map { meta, cram -> [ meta, cram, [] ] }
    | set { ch_cram_crai}

    fasta
    | map { meta, fasta -> fasta }
    | set { ch_fasta }

    SAMTOOLS_VIEW ( ch_cram_crai, ch_fasta, [] )
    ch_versions = ch_versions.mix ( SAMTOOLS_VIEW.out.versions.first() ) 


    // Calculate genome statistics
    FASTAWINDOWS ( fasta )
    ch_versions = ch_versions.mix ( FASTAWINDOWS.out.versions.first() )


    // Create genome windows file in BED format
    CREATE_BED ( FASTAWINDOWS.out.mononuc )
    ch_versions = ch_versions.mix ( CREATE_BED.out.versions.first() )

    
    // Calculate coverage
    SAMTOOLS_VIEW.out.bam
    | join ( SAMTOOLS_VIEW.out.csi )
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


    emit:
    freq     = FASTAWINDOWS.out.freq       // channel: [ val(meta), path(freq) ]
    mononuc  = FASTAWINDOWS.out.mononuc    // channel: [ val(meta), path(mononuc) ]
    bed      = CREATE_BED.out.bed          // channel: [ val(meta), path(bed) ]
    cov      = ch_coverage                 // channel: [ val(meta), path(regions.bed.gz) ]
    versions = ch_versions                 // channel: [ versions.yml ]
}
