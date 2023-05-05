//
// Generate static plots and summary from BlobDir
//

include { BLOBTOOLKIT_IMAGES  } from '../../modules/local/blobtoolkit/images'

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
    
    BLOBTOOLKIT_IMAGES ( blobdir, plots )
    ch_versions = ch_versions.mix( BLOBTOOLKIT_IMAGES.out.versions )
    

    emit:
    images   = BLOBTOOLKIT_IMAGES.out.png    // channel: [ val(meta), path(png)  ]
    versions = ch_versions                   // channel: [ versions.yml ]
}
