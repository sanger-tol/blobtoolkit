//
// Generate summary and static plots from blobdir
//

include { BLOBTOOLKIT_SUMMARY } from '../../modules/local/blobtoolkit/summary'
include { BLOBTK_IMAGES       } from '../../modules/local/blobtk/images'

workflow VIEW {
    take:
    blobdir     // channel: [ val(meta), path(blobdir) ]


    main:
    ch_versions = Channel.empty()


    //
    // Generate summary file
    //
    BLOBTOOLKIT_SUMMARY ( blobdir )
    ch_versions = ch_versions.mix ( BLOBTOOLKIT_SUMMARY.out.versions.first() )


    //
    // Generate static plots in png format
    //
    plots = [ "blob", "cumulative", "snail" ]

    BLOBTK_IMAGES ( blobdir, plots )
    ch_versions = ch_versions.mix( BLOBTK_IMAGES.out.versions )


    emit:
    summary  = BLOBTOOLKIT_SUMMARY.out.json  // channel: [ val(meta), path(json) ]
    images   = BLOBTK_IMAGES.out.png         // channel: [ val(meta), path(png) ]
    versions = ch_versions                   // channel: [ versions.yml ]
}
