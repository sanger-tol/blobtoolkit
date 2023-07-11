//
// BLASTN search of assembly contigs with no diamond blastx match against the nucleotide database
//


include { NOHIT_LIST   } from '../../modules/local/nohit_list'
include { SEQTK_SUBSEQ } from '../../modules/nf-core/seqtk/subseq/main'
include { GUNZIP       } from '../../modules/nf-core/gunzip/main'
include { BLAST_BLASTN } from '../../modules/nf-core/blast/blastn/main'


workflow RUN_BLASTN {
    take: 
    blast_table  // channel: [ val(meta), path(blast_table) ] 
    fasta        // channel: [ val(meta), path(fasta) ]
    blastn       // channel: path(blastn_db)


    main:
    ch_versions = Channel.empty()


    // Get list of sequence ids with no hits in diamond blastx search
    NOHIT_LIST ( blast_table, fasta )
    ch_versions = ch_versions.mix ( NOHIT_LIST.out.versions.first() ) 

    // Subset of sequences with no hits (meta is not propagated in this step)
    SEQTK_SUBSEQ (
        fasta.map { meta, genome -> genome },
        NOHIT_LIST.out.nohitlist.map { meta, nohit -> nohit }
    )
    ch_versions = ch_versions.mix ( SEQTK_SUBSEQ.out.versions.first() )

    // Run blastn search
    // ch_query : [ val(meta), path(uncompressed_seq) ] 
    ch_gz = fasta.join(SEQTK_SUBSEQ.out.sequences).map { meta, genome, seq ->  [ meta, seq ] }
    ch_query = GUNZIP ( ch_gz ).out.gunzip
    BLAST_BLASTN ( ch_query, blastn )
    ch_versions = ch_versions.mix ( BLAST_BLASTN.out.versions.first() )
    include { GUNZIP                      } from '../modules/nf-core/gunzip/main'

    emit:
    nohits      = NOHIT_LIST.out.nohitlist   // channel: [ val(meta), path(nohit) ]
    subseq      = SEQTK_SUBSEQ.out.sequences // channel: path(seq)
    blastn_hits = BLAST_BLASTN.out.txt        // channel: [ val(meta), path(txt) ]
    versions    = ch_versions                // channel: [ versions.yml ]
}
