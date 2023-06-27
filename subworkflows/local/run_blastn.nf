//
// BLASTN search of assembly contigs with no diamond blastx match against the nucleotide database
//


include { NOHIT_LIST   } from '../../modules/local/nohit_list'
include { SEQTK_SUBSEQ } from '../../modules/nf-core/seqtk/subseq/main'


workflow RUN_BLASTN {
    take: 
    blast_table    // channel: [ val(meta), path(blast_table) ] 
    fasta          // channel: [ val(meta), path(fasta) ]


    main:
    ch_versions = Channel.empty()


    // Get list of sequence ids with no hits in diamond blastx search
    NOHIT_LIST ( blast_table, fasta )
    ch_versions = ch_versions.mix ( NOHIT_LIST.out.versions.first() ) 

    // Subset of sequences with no hits
    SEQTK_SUBSEQ (
        fasta.map { meta, genome -> genome },
        NOHIT_LIST.out.nohitlist.map { meta, nohit -> nohit }
    )

    emit:
    nohits    = NOHIT_LIST.out.nohitlist   // channel: [ val(meta), path(freq) ]
    subseq    = SEQTK_SUBSEQ.out.sequences // channel: path() 
    versions  = ch_versions                // channel: [ versions.yml ]
}
