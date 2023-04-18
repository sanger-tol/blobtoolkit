//
//
// runs diamond_blastx taking BUSCO chunks as input
//
//

nextflow.enable.dsl = 2

include { CHUNK_FASTA_BY_BUSCO } from '../../modules/local/chunk_fasta_by_busco'
include { UNCHUNK_BLASTX       } from '../../modules/local/unchunk_blastx'
include { DIAMOND_BLASTX       } from '../../modules/nf-core/diamond/blastx/main'

workflow BUSCO_DIAMOND {
    take:

    // genome fasta
    fasta
    // first BUSCO lineage table
    busco_table
    // bed file
    bed

    main:

    ch_versions = Channel.empty()

    //
    // Chunk input fasta file
    //
    CHUNK_FASTA_BUSCO (
    fasta,
    busco_table,
    bed
    )
    ch_versions = ch_versions.mix(CHUNK_FASTA_BUSCO.out.versions)
    
    //
    // runs DIAMOND_BLASTX if fasta file from EXTRACT_BUSCO_GENE is not empty
    //
    
    // path to diamond 
    blastx_db = Channel.fromPath(params.diamondblastx_db) 
    // checks if output fasta is empty
    empty_fasta = CHUNK_FASTA_BY_BUSCO.out.chunks.map { meta,p -> file("$p").isEmpty() }    
    if ( empty_fasta == 'true' ) {
    dmd_input_ch = Channel.empty()
    }
    else {
    dmd_input_ch = CHUNK_FASTA_BY_BUSCO.out.chunks.combine(blastx_db)
    }
    // run diamond_blastx
    DIAMOND_BLASTX (
    dmd_input_ch.map { [it[0], it[1]] },
    dmd_input_ch.map { it[2] },
    "${params.blastx_outext}",
    "${params.blastx_cols}"
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
