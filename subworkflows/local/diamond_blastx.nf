//
// Runs diamond blastx
//

include { CHUNK_FASTA_BUSCO } from '../../modules/local/chunk_fasta_by_busco'
include { UNCHUNK_BLASTX    } from '../../modules/local/unchunk_blastx'
include { DIAMOND_BLASTX    } from '../../modules/nf-core/modules/diamond/blastx/main'

workflow CHUNK_BLASTX {
    take:

    // CHUNK_FASTA_BUSCO input

    /// File: Input fasta
    fasta
    /// File: Table for the first BUSCO lineage, it is an output from busco_diamond_blastp
    busco_table


    main:

    ch_versions = Channel.empty()

    //
    // Chunk input fasta file
    //
    CHUNK_FASTA_BUSCO (
    fasta,
    busco_table
    )
    ch_versions = ch_versions.mix(CHUNK_FASTA_BUSCO.out.versions.first())

    //
    // Runs diamond_blastx taking fasta chunks as input
    //
    DIAMOND_BLASTX (
    CHUNK_FASTA_BUSCO.out.chunks
    )
    ch_versions = ch_versions.mix(DIAMOND_BLASTX.out.versions)

    //
    // Unchunk_blastx results
    //
    UNCHUNK_BLASTX (
    DIAMOND_BLASTX.out.blast
    )
    ch_versions = ch_versions.mix(UNCHUNK_BLASTX.out.versions)

    emit:

    //  UNCHUNK_BLASTX output
    reference_proteomes  = UNCHUNK_BLASTX.out.proteomes
    // tool versions
    versions = ch_versions
}
