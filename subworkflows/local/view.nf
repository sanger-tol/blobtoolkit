//
// Create summary and plots from blobdir
//

include { BLOBTOOLKIT_SUMMARY } from '../../modules/local/blobtoolkit/summary'

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


    emit:
    summary = BLOBTOOLKIT_SUMMARY.out.json   // channel: [ val(meta), path(json) ]
    versions = ch_versions                   // channel: [ versions.yml ]
}
