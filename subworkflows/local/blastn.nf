//
// Runs blast/blastn
//

include { EXTRACT_NOHIT_FASTA                    } from '../../modules/local/extract_nohit_fasta'
include { CHUNK_FASTA_BUSCO as CHUNK_NOHIT_FASTA } from '../../modules/local/chunk_fasta_by_busco'
include { RUN_BLASTN                             } from '../../modules/local/run_blastn'
include { UNCHUNK_BLASTN                         } from '../../modules/local/unchunk_blastn'

workflow BLASTN {
    take:

    // EXTRACT_NOHIT_FASTA input
    // File: genome fasta File
    fasta
    // File: blastx results from CHUNK_BLASTX (diamond_blastx) subworkflow
    blastx_table
    // Value: evalue to filter input fasta file from blastx_table
    blastx_evalue

    // CHUNK_NOHIT_FASTA input
    /// File: Input fasta from EXTRACT_NOHIT_FASTA
    /// File: Table for the first BUSCO lineage, it is an output from busco_diamond_blastp
    busco_table
    /// Value: minimum chunk size for splitting long sequences, default = 100000
    blast_chunk
    /// Value: overlap length for splitting long sequences, default = 0
    blast_chunk_overlap
    /// Value: minimum chunk size for splitting long sequences, default = 10
    blast_max_chunks
    /// Value: minimum chunk size for splitting long sequences, default = 1000
    blast_min_length

    // RUN_BLASTN
    //  File: Input fasta file is the output from CHUNK_NOHIT_FASTA
    /// Path: database for BLASTN:
    blastn_db
    /// Value: evalue for BLASTN
    blastn_evalue


    // UNCHUNK_BLASTN input
    //  File: Input file is the output from RUN_BLASTN
    /// Value: max_target_seqs, same as in RUN_BLASTN
    max_target_seqs


    main:

    ch_versions = Channel.empty()
    // not uses as it is declared within each module
    // name = fasta.simpleName()

    //
    // Extract sequences with no blastx hits into a separate file
    //
    EXTRACT_NOHIT_FASTA (
    fasta,
    blastx_table,
    blastx_evalue
    )
    ch_versions = ch_versions.mix(EXTRACT_NOHIT_FASTA.out.versions.first())

    //
    // Chunk no-hit fasta file
    //
    CHUNK_NOHIT_FASTA (
    EXTRACT_NOHIT_FASTA.out.nohit_fasta,
    busco_table,
    blast_chunk,
    blast_chunk_overlap,
    blast_max_chunks,
    blast_min_length
    )
    ch_versions = ch_versions.mix(CHUNK_NOHIT_FASTA.out.versions)

    //
    // Runs blastn
    //
    RUN_BLASTN (
    CHUNK_NOHIT_FASTA.out.chunks,
    blastn_db,
    blastn_evalue,
    max_target_seqs
    )
    ch_versions = ch_versions.mix(RUN_BLASTN.out.versions)

    //
    // Unchunk_blastn results
    //
    UNCHUNK_BLASTN (
    RUN_BLASTN.out.blastn_out,
    max_target_seqs
    )
    ch_versions = ch_versions.mix(UNCHUNK_BLASTN.out.versions)

    emit:

    //  UNCHUNK_BLASTN output
    blastn_unchunk  = UNCHUNK_BLASTN.out.blastn_out
    // tool versions
    versions = ch_versions
}
