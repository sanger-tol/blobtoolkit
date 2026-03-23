//
// Prepare genome for downstream processing
//

include { GUNZIP                } from '../../modules/nf-core/gunzip/main'

include { REPEAT_MASKING        } from '../sanger-tol/repeat_masking/main'


workflow PREPARE_GENOME {
    take:
    genome    // channel: [ meta, path(fasta) ]


    main:
    ch_versions = channel.empty()

    //
    // LOGIC: Identify the compressed files
    //
    ch_genomes_for_gunzip = genome
        .branch { _meta, fasta ->
            gunzip: fasta.name.endsWith( ".gz" )
            skip: true
        }


    //
    // MODULE: Decompress compressed FASTA files
    //
    GUNZIP ( ch_genomes_for_gunzip.gunzip )


    //
    // LOGIC: Extract the genome size for decision making downstream
    //
    ch_genome = ch_genomes_for_gunzip.skip
        .mix(GUNZIP.out.gunzip)
        .map { meta, fa -> [ meta + [genome_size: fa.size()], fa] }


    //
    // SUBWORKFLOW: Mask the genome if needed
    //
    if ( params.mask ) {
        REPEAT_MASKING ( ch_genome )

        ch_fasta = REPEAT_MASKING.out.repeat_intervals
    } else {
        ch_fasta = ch_genome
    }


    emit:
    genome   = ch_fasta            // channel: [ meta, path(genome) ]
    versions = ch_versions         // channel: [ versions.yml ]
}
