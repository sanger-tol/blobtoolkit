//
// Create BlobTools dataset
//

include { BLOBTOOLKIT_CHUNK_BUSCO } from '../../modules/local/blobtoolkit/chunk_busco'
include { DIAMOND_BLASTX          } from '../../modules/nf-core/diamond/blastx/main'

workflow RUN_DIAMOND_BLASTX {
    take:
    fasta      // channel: [ val(meta), path(fasta) ]
    table      // channel: [ val(meta), path(busco_table) ]
    blastx     // channel: path(blastx_db)
    outext     // channel: val(out_format)
    cols       // channel: val(column_names)


    main:
    ch_versions = Channel.empty()


    //
    // Create metadata summary file
    //
    BLOBTOOLKIT_CHUNK_BUSCO ( fasta, table )
    ch_versions = ch_versions.mix ( BLOBTOOLKIT_CHUNK_BUSCO.out.versions.first() )

    //
    // Run diamond_blastx
    //
    DIAMOND_BLASTX ( BLOBTOOLKIT_CHUNK_BUSCO.out.chunks, blastx, outext, cols)


    emit:
    fasta_chunks = BLOBTOOLKIT_CHUNK_BUSCO.out.chunks  // channel: [ val(meta), path(chunks) ]
    versions     = ch_versions                         // channel: [ versions.yml ]
}
