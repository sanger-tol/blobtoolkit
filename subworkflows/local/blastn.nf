//
// Runs blast/blastn
//

include { GET_NOHIT_LIST                         } from '../../modules/local/get_nohit_list'
include { EXTRACT_NOHIT_FASTA                    } from '../../modules/local/extract_nohit_fasta'
include { CHUNK_FASTA_BUSCO as CHUNK_NOHIT_FASTA } from '../../modules/local/chunk_fasta_by_busco'
include { RUN_BLASTN                             } from '../../modules/local/run_blastn'
include { UNCHUNK_BLASTN                         } from '../../modules/local/unchunk_blastn'

workflow BLASTN {
    take:

    // File: genome fasta File
    fasta

    // File: blastx results from CHUNK_BLASTX (diamond_blastx) subworkflow
    blastx_table

    // CHUNK_NOHIT_FASTA input
    busco_table

    /// Path: database for BLASTN:
    blastn_db


    main:

    ch_versions = Channel.empty()

    //
    // Extract sequences with no blastx hits into a separate file
    //

    GET_NOHIT_LIST (
    fasta,
    blastx,
    params.evalue
    )
    ch_versions = ch_versions.mix(GET_NOHIT_LIST.out.versions)

    EXTRACT_NOHIT_FASTA (
    fasta,
    GET_NOHIT_LIST.out.nohit_list
    )
    ch_versions = ch_versions.mix(EXTRACT_NOHIT_FASTA.out.versions)

    //
    // Chunk no-hit fasta file
    //
    CHUNK_NOHIT_FASTA (
    fasta,
    busco_table
    )
    ch_versions = ch_versions.mix(CHUNK_NOHIT_FASTA.out.versions)

    //
    // Runs blastn
    //
    RUN_BLASTN (
    CHUNK_NOHIT_FASTA.out.chunks,
    )
    ch_versions = ch_versions.mix(RUN_BLASTN.out.versions)

    //
    // Unchunk_blastn results
    //
    UNCHUNK_BLASTN (
    RUN_BLASTN.out.blastn_out
    )
    ch_versions = ch_versions.mix(UNCHUNK_BLASTN.out.versions)

    emit:

    //  UNCHUNK_BLASTN output
    blastn_unchunk  = UNCHUNK_BLASTN.out.blastn_out
    // tool versions
    versions = ch_versions
}
