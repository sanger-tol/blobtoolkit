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
    //
    // MODULE: CREATE BLOBTOOLS DATASET FILES
    //
    BLOBTOOLKIT_CREATEBLOBDIR (
        windowstats,
        busco,
        blastp,
        config,
        taxdump
    )


    //
    // MODULE: UPDATE BLOBTOOLS DATASET FILES
    //
    BLOBTOOLKIT_UPDATEBLOBDIR (
        BLOBTOOLKIT_CREATEBLOBDIR.out.blobdir,
        syn_tsv,
        cat_tsv,
        blastx,
        blastn,
        taxdump
    )


    emit:
    blobdir  = BLOBTOOLKIT_UPDATEBLOBDIR.out.blobdir  // channel: [ val(meta), path(dir) ]
}
