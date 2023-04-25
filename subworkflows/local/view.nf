//
// Generate static plots and summary from BlobDir
//

include { BLOBTOOLKIT_IMAGES  } from '../../modules/local/blobtoolkit/images'
include { BLOBTOOLKIT_SUMMARY } from '../../modules/local/blobtoolkit/summary'

workflow VIEW {
    take:
    blobdir // channel: [ val(meta), path(blobdir) ]

    main:
    ch_versions = Channel.empty()

    //
    // Generate static plots in png format
    //
    
    plots = [
    "--view blob --param plotShape=circle",
    "--view blob --param plotShape=hex",
    "--view blob --param plotShape=square",
    "--view blob --param plotShape=kite",
    "--view cumulative",
    "--view snail"
    ]
    
    GENERATE_IMAGES ( blobdir, plots )
    ch_versions = ch_versions.mix( GENERATE_IMAGES.out.versions )
    

    //  
    // Generate summary
    //
    GENERATE_SUMMARY ( blobdir )
    ch_versions = ch_versions.mix( GENERATE_SUMMARY.out.versions )


    emit:
    images   = GENERATE_IMAGES.out.png    // channel: [ val(meta), path(png)  ]
    summary  = GENERATE_SUMMARY.out.json  // channel: [ val(meta), path(json) ]
    versions = ch_versions                // channel: [ versions.yml ]
}
