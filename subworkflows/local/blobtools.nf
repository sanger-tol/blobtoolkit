//
// Create BlobTools dataset
//

include { BLOBTOOLKIT_METADATA } from '../../modules/local/blobtoolkit/metadata'
include { CREATE_BLOBDIR       } from '../../modules/local/create_blobdir'

workflow BLOBTOOLS {
    take:
    config      // channel: [ val(meta), path(config) ]
    windowstats // channel: [ val(meta), path(window_stats_tsvs) ]
    busco       // channel: [ val(meta), path(full_table) ]
    blastp      // channel: [ val(meta), path(txt) ]
    taxdump     // channel: path(taxdump_db)


    main:
    ch_versions = Channel.empty()


    //
    // Create metadata summary file
    //
    BLOBTOOLKIT_METADATA ( config )
    ch_versions = ch_versions.mix ( BLOBTOOLKIT_METADATA.out.versions.first() )


    //  
    // Create Blobtools dataset files
    //
    CREATE_BLOBDIR ( windowstats, busco, blastp, BLOBTOOLKIT_METADATA.out.yaml, taxdump )
    ch_versions = ch_versions.mix ( CREATE_BLOBDIR.out.versions.first() )


    emit:
    metadata = BLOBTOOLKIT_METADATA.out.yaml // channel: [ val(meta), path(yaml) ]
    blobdir  = CREATE_BLOBDIR.out.blobdir       // channel: [ val(meta), path(dir) ]
    versions = ch_versions                      // channel: [ versions.yml ]
}
