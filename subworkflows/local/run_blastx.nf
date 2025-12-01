//
// Diamond blastx search of assembly contigs against the UniProt reference proteomes
//

include { BLOBTOOLKIT_CHUNK   } from '../../modules/local/blobtoolkit/chunk'
include { BLOBTOOLKIT_UNCHUNK } from '../../modules/local/blobtoolkit/unchunk'
include { DIAMOND_BLASTX      } from '../../modules/nf-core/diamond/blastx/main'

workflow RUN_BLASTX {
    take:
    fasta      // channel: [ val(meta), path(fasta) ]
    table      // channel: [ val(meta), path(busco_table) ]
    blastx     // channel: [ val(meta), path(blastx_db) ]
    taxon_id   // channel: val(taxon_id)


    main:
    ch_versions = Channel.empty()


    //
    // Split the sequences
    //
    BLOBTOOLKIT_CHUNK ( fasta, table )
    ch_versions = ch_versions.mix ( BLOBTOOLKIT_CHUNK.out.versions.first() )


    //
    // Run diamond_blastx
    //
    // Hardocded to match the format expected by blobtools
    def outext = 'txt'
    def cols   = 'qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
    DIAMOND_BLASTX ( BLOBTOOLKIT_CHUNK.out.chunks, blastx, outext, cols, taxon_id )
    ch_versions = ch_versions.mix ( DIAMOND_BLASTX.out.versions.first() )


    //
    // Unchunk chunked blastx results
    //
    BLOBTOOLKIT_UNCHUNK ( DIAMOND_BLASTX.out.txt )
    ch_versions = ch_versions.mix ( BLOBTOOLKIT_UNCHUNK.out.versions.first() )


    emit:
    blastx_out = BLOBTOOLKIT_UNCHUNK.out.blast_out  // channel: [ val(meta), path(blastx_out) ]
    versions   = ch_versions                        // channel: [ versions.yml ]
}
