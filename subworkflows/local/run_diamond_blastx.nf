//
// Create BlobTools dataset
//

include { BLOBTOOLKIT_CHUNK_BUSCO } from '../../modules/local/blobtoolkit/chunk_busco'

workflow RUN_DIAMOND_BLASTX {
    take:
    fasta      // channel: [ val(meta), path(fasta) ]
    table      // channel: [ val(meta), path(busco_table) ]


    main:
    ch_versions = Channel.empty()


    //
    // Create metadata summary file
    //
    BLOBTOOLKIT_CHUNK_BUSCO ( fasta, table )
    ch_versions = ch_versions.mix ( BLOBTOOLKIT_CHUNK_BUSCO.out.versions.first() ) 


    emit:
    fasta_chunks = BLOBTOOLKIT_CHUNK_BUSCO.out.chunks  // channel: [ val(meta), path(chunks) ]
    versions     = ch_versions                         // channel: [ versions.yml ]
}
