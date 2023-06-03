//
// Create BlobTools dataset
//

include { BLOBTOOLKIT_CHUNK_BUSCO } from '../../modules/local/blobtoolkit/chunk_busco'
include { BLOBTOOLKIT_UNCHUNK     } from '../../modules/local/blobtoolkit/unchunk_blastx'
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

    //
    // Unchunk chunked blastx results
    //
    BLOBTOOLKIT_UNCHUNK ( DIAMOND_BLASTX.out.txt )


    emit:
    blastx_out = BLOBTOOLKIT_UNCHUNK.out.blastx  // channel: [ val(meta), path(blastx_out) ]
    versions   = ch_versions                     // channel: [ versions.yml ]
}
