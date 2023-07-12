//
// BLASTN search of assembly contigs with no diamond blastx match against the nucleotide database
//


include { NOHIT_LIST          } from '../../modules/local/nohit_list'
include { SEQTK_SUBSEQ        } from '../../modules/nf-core/seqtk/subseq/main'
include { GUNZIP              } from '../../modules/nf-core/gunzip/main'
include { BLOBTOOLKIT_CHUNK   } from '../../modules/local/blobtoolkit/chunk'
include { BLAST_BLASTN        } from '../../modules/nf-core/blast/blastn/main'


workflow RUN_BLASTN {
    take: 
    blast_table  // channel: [ val(meta), path(blast_table) ] 
    fasta        // channel: [ val(meta), path(fasta) ]
    blastn       // channel: path(blastn_db)


    main:
    ch_versions = Channel.empty()


    // Extract no hits fasta
    // Get list of sequence ids with no hits in diamond blastx search
    NOHIT_LIST ( blast_table, fasta )
    ch_versions = ch_versions.mix ( NOHIT_LIST.out.versions.first() ) 
    // Subset of sequences with no hits (meta is not propagated in this step)
    SEQTK_SUBSEQ (
        fasta.map { meta, genome -> genome },
        NOHIT_LIST.out.nohitlist.map { meta, nohit -> nohit }
    )
    ch_versions = ch_versions.mix ( SEQTK_SUBSEQ.out.versions.first() )
    
    
    //  Split long contigs into chunks 
    // add meta to fasta subset channel: [ val(meta), path(compressed_fasta) ]
    ch_gz = fasta.combine(SEQTK_SUBSEQ.out.sequences).map { meta, genome, seq ->  [ meta, seq ] }
    // uncompress fasta
    GUNZIP ( ch_gz )
    // create chunks
    BLOBTOOLKIT_CHUNK ( GUNZIP.out.gunzip, [[],[]] )
    ch_versions = ch_versions.mix ( BLOBTOOLKIT_CHUNK.out.versions.first() )


    // Run blastn search
    // channel to query fasta: [ val(meta), path(uncompressed_fasta) ] 
    BLAST_BLASTN ( BLOBTOOLKIT_CHUNK.out.chunks, blastn )
    ch_versions = ch_versions.mix ( BLAST_BLASTN.out.versions.first() )


    emit:
    nohits      = NOHIT_LIST.out.nohitlist   // channel: [ val(meta), path(nohit) ]
    subseq      = SEQTK_SUBSEQ.out.sequences // channel: path(seq)
    blastn_hits = BLAST_BLASTN.out.txt       // channel: [ val(meta), path(txt) ]
    versions    = ch_versions                // channel: [ versions.yml ]
}
