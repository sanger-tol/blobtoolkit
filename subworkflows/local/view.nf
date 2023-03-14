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
    GENERATE_IMAGES (
      blobdir
    )
    
    //
    // Generate summary
    //
    // channel: [meta, blobdir, png]
    summary_input = blobdir.combine(GENERATE_IMAGES.out.png, by:0)
    GENERATE_SUMMARY (
      summary_input.map { [it[0],it[1]] },
      summary_input.map { [it[0],it[2]] }
    )
    
    emit:

    // png file
    png_img = GENERATE_IMAGES.out.png

    // summary json file 
    summary_json = GENERATE_SUMMARY.out.json

    // tool versions
    versions = ch_versions
}
