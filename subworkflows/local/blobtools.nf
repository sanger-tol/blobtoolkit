//
// Create BlobTools dataset
//

include { BLOBTOOLKIT_CREATEBLOBDIR } from '../../modules/local/blobtoolkit/createblobdir'
include { BLOBTOOLKIT_UPDATEBLOBDIR } from '../../modules/local/blobtoolkit/updateblobdir'

workflow BLOBTOOLS {
    take:
    config      // channel: [ val(meta), path(config) ]
    syn_tsv     // channel: [ val(meta), [path(tsv)] ]
    cat_tsv     // channel: [ val(meta), [path(tsv)] ]
    windowstats // channel: [ val(meta), path(window_stats_tsvs) ]
    busco       // channel: [ val(meta), path(full_table) ]
    blastp      // channel: [ val(meta), path(txt) ]
    blastx      // channel: [ val(meta), path(txt) ]
    blastn      // channel: [ val(meta), path(txt) ]
    taxdump     // channel: path(taxdump_db)


    main:
    ch_versions = Channel.empty()


    //
    // Create Blobtools dataset files
    //
    BLOBTOOLKIT_CREATEBLOBDIR ( windowstats, busco, blastp, config, taxdump )
    ch_versions = ch_versions.mix ( BLOBTOOLKIT_CREATEBLOBDIR.out.versions.first() )


    //
    // Update Blobtools dataset files
    //
    BLOBTOOLKIT_UPDATEBLOBDIR ( BLOBTOOLKIT_CREATEBLOBDIR.out.blobdir, syn_tsv, cat_tsv, blastx, blastn, taxdump )
    ch_versions = ch_versions.mix ( BLOBTOOLKIT_UPDATEBLOBDIR.out.versions.first() )


    emit:
    blobdir  = BLOBTOOLKIT_UPDATEBLOBDIR.out.blobdir  // channel: [ val(meta), path(dir) ]
    versions = ch_versions                            // channel: [ versions.yml ]
}
