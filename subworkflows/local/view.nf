//
//
// Generate static plots and summary from BlobDir
//
//

nextflow.enable.dsl = 2

include { GENERATE_IMAGES  } from '../../modules/local/generate_images'
include { GENERATE_SUMMARY } from '../../modules/local/generate_summary'

workflow VIEW {
    take:

    //  Tuple [meta, blobdir]:
    blobdir

    main:

    ch_versions = Channel.empty()
    
    //
    // Generate static plot (png format)
    //
    
    // list to select plots
    plots = [
      "--view blob --param plotShape=circle",
      "--view blob --param plotShape=hex",
      "--view blob --param plotShape=square",
      "--view blob --param plotShape=kite",
      "--view cumulative",
      "--view snail"
    ]
    
    GENERATE_IMAGES (
      blobdir,
      ["--view cumulative"]
    )
    ch_versions = ch_versions.mix(GENERATE_IMAGES.out.versions)
    
    //
    // Generate summary
    //
    // channel: [meta, blobdir, png]
    GENERATE_SUMMARY (
      blobdir
    )
    ch_versions = ch_versions.mix(GENERATE_SUMMARY.out.versions)

    emit:

    // png file
    png_img = GENERATE_IMAGES.out.png

    // summary json file 
    summary_json = GENERATE_SUMMARY.out.json

    // tool versions
    versions = ch_versions
}
