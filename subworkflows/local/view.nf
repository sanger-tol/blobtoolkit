//
// Generate summary and static plots from blobdir
//

include { BLOBTOOLKIT_SUMMARY } from '../../modules/local/blobtoolkit/summary'
include { BLOBTK_IMAGES       } from '../../modules/local/blobtk/images'

workflow VIEW {
    take:
    blobdir     // channel: [ val(meta), path(blobdir) ]


    main:
    ch_versions = channel.empty()


    //
    // Generate summary file
    //
    BLOBTOOLKIT_SUMMARY ( blobdir )
    ch_versions = ch_versions.mix ( BLOBTOOLKIT_SUMMARY.out.versions.first() )


    //
    // Generate static plots in png/svg format
    //
    plots = [ "blob", "cumulative", "snail" ]

    BLOBTK_IMAGES ( blobdir, plots, params.image_format )
    ch_versions = ch_versions.mix( BLOBTK_IMAGES.out.versions )

    ch_images = BLOBTK_IMAGES.out.png.mix(BLOBTK_IMAGES.out.svg)

    emit:
    summary  = BLOBTOOLKIT_SUMMARY.out.json  // channel: [ val(meta), path(json) ]
    images   = ch_images                     // channel: [ val(meta), path(png/svg) ]
    versions = ch_versions                   // channel: [ versions.yml ]
}
