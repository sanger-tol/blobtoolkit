//
// Final edits to the blobdir
//

include { BLOBTOOLKIT_UPDATEMETA } from '../../modules/local/blobtoolkit/updatemeta'
include { COMPRESSBLOBDIR        } from '../../modules/local/compressblobdir'

workflow FINALISE_BLOBDIR {
    take:
    blobdir     // channel: [ val(meta), path(blobdir) ]
    software    // channel: [ val(meta), path(software_yml) ]
    summary     // channel: [ val(meta), path(summary_json) ]


    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Update the software listed in the meta json file
    //
    BLOBTOOLKIT_UPDATEMETA ( blobdir, software )

    //
    // MODULE: Compress all the json files
    //
    COMPRESSBLOBDIR ( blobdir, summary, BLOBTOOLKIT_UPDATEMETA.out.json )
    ch_versions = ch_versions.mix ( COMPRESSBLOBDIR.out.versions.first() )


    emit:
    blobdir  = COMPRESSBLOBDIR.out.blobdir  // channel: [ val(meta), path(json) ]
    versions = ch_versions                              // channel: [ versions.yml ]
}
