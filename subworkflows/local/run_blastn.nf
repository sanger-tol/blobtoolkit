//
// BLASTN search of assembly contigs with no diamond blastx match against the nucleotide database
//


include { NOHIT_LIST             } from '../../modules/local/nohit_list'
include { SEQTK_SUBSEQ           } from '../../modules/nf-core/seqtk/subseq/main'
include { GUNZIP                 } from '../../modules/nf-core/gunzip/main'
include { BLOBTOOLKIT_CHUNK      } from '../../modules/local/blobtoolkit/chunk'
include { BLASTN as BLASTN_TAXON } from '../../modules/local/blastn'
include { BLASTN                 } from '../../modules/local/blastn'
include { BLOBTOOLKIT_UNCHUNK    } from '../../modules/local/blobtoolkit/unchunk'


workflow RUN_BLASTN {
    take: 
    blast_table  // channel: [ val(meta), path(blast_table) ] 
    fasta        // channel: [ val(meta), path(fasta) ]
    blastn       // channel: path(blastn_db)
    taxon_id     // channel: val(taxon_id)


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
    // run blastn excluding taxon_id
    BLASTN_TAXON ( BLOBTOOLKIT_CHUNK.out.chunks, blastn, taxon_id )
    // check if blastn output table is empty
    BLASTN_TAXON.out.txt
    | map { meta, txt -> txt.isEmpty() }
    | set { is_txt_empty }
    // repeat the blastn search without excluding taxon_id
    if ( is_txt_empty ) {
    BLASTN ( BLOBTOOLKIT_CHUNK.out.chunks, blastn, [] )
    ch_blastn_txt = BLASTN.out.txt
    }
    else {
    ch_blastn_txt = BLASTN_TAXON.out.txt
    }

    ch_versions = ch_versions.mix ( BLASTN.out.versions.first() )


    // Unchunk chunked blastn results
    BLOBTOOLKIT_UNCHUNK ( ch_blastn_txt )
    ch_versions = ch_versions.mix ( BLOBTOOLKIT_UNCHUNK.out.versions.first() )


    emit:
    blastn_out = BLOBTOOLKIT_UNCHUNK.out.blast_out  // channel: [ val(meta), path(blastn_out) ]
    versions   = ch_versions                        // channel: [ versions.yml ]
}
