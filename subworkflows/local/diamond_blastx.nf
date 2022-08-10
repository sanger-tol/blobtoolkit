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
    /// Value: minimum chunk size for splitting long sequences, default = 100000
    blast_chunk
    /// Value: overlap length for splitting long sequences, default = 0
    blast_chunk_overlap
    /// Value: minimum chunk size for splitting long sequences, default = 10
    blast_max_chunks
    /// Value: minimum chunk size for splitting long sequences, default = 1000
    blast_min_length

    // DIAMOND_BLASTX
    //  File: Input fasta file is the output from CHUNK_FASTA_BUSCO
    //  args: evalue, max_target_seqs
    /// Path: Directory containing the diamond database (for regions):
    diamonddb_blastx
    /// Value: Specify the type of output file to be generated, `txt` corresponds to to BLAST tabular format:
    outext
    /// Value: pace separated list of DIAMOND tabular BLAST output keywords:
    /// "qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    blastx_cols

    // UNCHUNK_BLASTX input
    //  File: Input file is the raw proteome (blast?) from DIAMOND_BLASTX
    /// Value: max_target_seqs, same as in DIAMOND_BLASTX args
    max_target_seqs


    main:

    ch_versions = Channel.empty()
    name = fasta.simpleName()

    //
    // Chunk input fasta file
    //
    CHUNK_FASTA_BUSCO (
    fasta,
    busco_table,
    blast_chunk,
    blast_chunk_overlap,
    blast_max_chunks,
    blast_min_length
    )
    ch_versions = ch_versions.mix(CHUNK_FASTA_BUSCO.out.versions.first())

    //
    // Runs diamond_blastx taking fasta chunks as input
    //
    DIAMOND_BLASTX (
    [ [ id:name ], CHUNK_FASTA_BUSCO.out.chunks ],
    diamonddb_blastx,
    outext,
    blastx_cols
    )
    ch_versions = ch_versions.mix(DIAMOND_BLASTX.out.versions)

    //
    // Unchunk_blastx results
    //
    UNCHUNK_BLASTX (
    DIAMOND_BLASTX.out.blast,
    max_target_seqs
    )
    ch_versions = ch_versions.mix(UNCHUNK_BLASTX.out.versions)

    emit:

    //  UNCHUNK_BLASTX outputs?
    reference_proteomes  = UNCHUNK_BLASTX.out.proteomes
    // tool versions
    versions = ch_versions
}
