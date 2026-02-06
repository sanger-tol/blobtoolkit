//
// BLASTN search of assembly contigs with no diamond blastx match against the nucleotide database
//


include { NOHIT_LIST                   } from '../../modules/local/nohit_list'
include { SEQTK_SUBSEQ                 } from '../../modules/nf-core/seqtk/subseq/main'
include { BLOBTOOLKIT_CHUNK            } from '../../modules/local/blobtoolkit/chunk'
include { BLAST_BLASTN as BLASTN_TAXON } from '../../modules/nf-core/blast/blastn/main'
include { BLAST_BLASTN                 } from '../../modules/nf-core/blast/blastn/main'
include { BLOBTOOLKIT_UNCHUNK          } from '../../modules/local/blobtoolkit/unchunk'


workflow RUN_BLASTN {
    take:
    blast_table  // channel: [ val(meta), path(blast_table) ]
    fasta        // channel: [ val(meta), path(fasta) ]
    blastn       // channel: [ val(meta), path(blastn_db) ]
    taxon_id     // channel: val(taxon_id)


    main:
    ch_versions = channel.empty()


    // Extract no hits fasta
    // Get list of sequence ids with no hits in diamond blastx search
    NOHIT_LIST ( blast_table, fasta )
    ch_versions = ch_versions.mix ( NOHIT_LIST.out.versions.first() )


    //
    // MODULE: Subset of sequences with no hits
    //
    SEQTK_SUBSEQ (
        fasta,
        NOHIT_LIST.out.nohitlist.map { _meta, nohit -> nohit } . filter { file -> file.size() > 0 }
    )
    ch_versions = ch_versions.mix ( SEQTK_SUBSEQ.out.versions.first() )


    //  Split long contigs into chunks
    // create chunks
    BLOBTOOLKIT_CHUNK ( SEQTK_SUBSEQ.out.sequences, [[],[]] )
    ch_versions = ch_versions.mix ( BLOBTOOLKIT_CHUNK.out.versions.first() )


    // Check that there are still sequences left after chunking (which excludes masked regions)
    BLOBTOOLKIT_CHUNK.out.chunks
    | filter { _meta, file -> file.size() > 0 }
    | set { ch_chunks }

    // Run blastn search
    if (params.skip_taxon_filtering) {

        // skip BLASTN_TAXON
        ch_blast_blastn_input = ch_chunks

        // fake ch_blastn_taxon_out.not_empty
        ch_blastn_taxon_out = [
            not_empty: channel.empty()
        ]

    } else {

        // run blastn excluding taxon_id
        BLASTN_TAXON ( ch_chunks, blastn, taxon_id )
        ch_versions = ch_versions.mix ( BLASTN_TAXON.out.versions.first() )

        // check if blastn output table is empty
        BLASTN_TAXON.out.txt
        | branch { _meta, txt ->
            empty:     txt.isEmpty()
            not_empty: true
        }
        | set { ch_blastn_taxon_out }

        // repeat the blastn search without excluding taxon_id
        ch_blastn_taxon_out.empty
        | join ( ch_chunks )
        | map { meta, _txt, file -> [meta, file] }
        | set { ch_blast_blastn_input }

    }

    BLAST_BLASTN ( ch_blast_blastn_input, blastn, [] )
    ch_versions = ch_versions.mix ( BLAST_BLASTN.out.versions.first() )

    BLAST_BLASTN.out.txt
    | mix( ch_blastn_taxon_out.not_empty )
    | set { ch_blastn_txt }


    //
    // MODULE: Unchunk chunked blastn results
    //
    BLOBTOOLKIT_UNCHUNK ( ch_blastn_txt )
    ch_versions = ch_versions.mix ( BLOBTOOLKIT_UNCHUNK.out.versions.first() )


    emit:
    blastn_out = BLOBTOOLKIT_UNCHUNK.out.blast_out  // channel: [ val(meta), path(blastn_out) ]
    versions   = ch_versions                        // channel: [ versions.yml ]
}
