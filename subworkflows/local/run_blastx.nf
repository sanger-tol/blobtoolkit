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
    //
    // Split the sequences
    //
    BLOBTOOLKIT_CHUNK ( fasta, table )


    //
    // Run diamond_blastx
    //
    // Hardocded to match the format expected by blobtools
    def outext = 'txt'
    def cols   = 'qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
    DIAMOND_BLASTX ( BLOBTOOLKIT_CHUNK.out.chunks, blastx, outext, cols, taxon_id )


    //
    // Unchunk chunked blastx results
    //
    BLOBTOOLKIT_UNCHUNK ( DIAMOND_BLASTX.out.txt )


    emit:
    blastx_out = BLOBTOOLKIT_UNCHUNK.out.blast_out  // channel: [ val(meta), path(blastx_out) ]
}
