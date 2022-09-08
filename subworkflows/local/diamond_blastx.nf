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


    main:

    ch_versions = Channel.empty()

    // this is the string used to name all intermediate and final output files
    name = fasta.map { f -> f.simpleName }

    //
    // Chunk input fasta file
    //
    CHUNK_FASTA_BUSCO (
    name,
    fasta,
    busco_table
    )
    ch_versions = ch_versions.mix(CHUNK_FASTA_BUSCO.out.versions.first())

    //
    // Runs diamond_blastx taking fasta chunks as input
    //
    DIAMOND_BLASTX (
    CHUNK_FASTA_BUSCO.out.chunks.map { fa -> [ [id: fa.baseName ], fa ] }, // Add meta data using the file's basename as id,
    diamonddb_blastx,
    outext,
    blastx_cols
    )
    ch_versions = ch_versions.mix(DIAMOND_BLASTX.out.versions)

    //
    // Unchunk_blastx results
    //
    UNCHUNK_BLASTX (
    DIAMOND_BLASTX.out.blast.map { f -> f.baseName },
    DIAMOND_BLASTX.out.blast
    )
    ch_versions = ch_versions.mix(UNCHUNK_BLASTX.out.versions)

    emit:

    //  UNCHUNK_BLASTX output
    reference_proteomes  = UNCHUNK_BLASTX.out.proteomes
    // tool versions
    versions = ch_versions
}
