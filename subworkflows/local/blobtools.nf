//
// Create BlobTools dataset
//

include { BLOBTOOLKIT_METADATA      } from '../../modules/local/blobtoolkit/metadata'
include { BLOBTOOLKIT_CREATEBLOBDIR } from '../../modules/local/blobtoolkit/createblobdir'
include { BLOBTOOLKIT_UPDATEBLOBDIR } from '../../modules/local/blobtoolkit/updateblobdir'

workflow BLOBTOOLS {
    take:
    config      // channel: [ val(meta), path(config) ]
    windowstats // channel: [ val(meta), path(window_stats_tsvs) ]
    busco       // channel: [ val(meta), path(full_table) ]
    blastp      // channel: [ val(meta), path(txt) ]
    blastx      // channel: [ val(meta), path(txt) ]
    blastn      // channel: [ val(meta), path(txt) ]
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
    BLOBTOOLKIT_CREATEBLOBDIR ( windowstats, busco, blastp, BLOBTOOLKIT_METADATA.out.yaml, taxdump )
    ch_versions = ch_versions.mix ( BLOBTOOLKIT_CREATEBLOBDIR.out.versions.first() )


    //
    // Update Blobtools dataset files
    //
    BLOBTOOLKIT_UPDATEBLOBDIR ( BLOBTOOLKIT_CREATEBLOBDIR.out.blobdir, blastx, blastn, taxdump )
    ch_versions = ch_versions.mix ( BLOBTOOLKIT_UPDATEBLOBDIR.out.versions.first() )


    emit:
    metadata = BLOBTOOLKIT_METADATA.out.yaml          // channel: [ val(meta), path(yaml) ]
    blobdir  = BLOBTOOLKIT_UPDATEBLOBDIR.out.blobdir  // channel: [ val(meta), path(dir) ]
    versions = ch_versions                            // channel: [ versions.yml ]
}
