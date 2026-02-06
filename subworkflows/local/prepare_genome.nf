//
// Prepare genome for downstream processing
//

include { GUNZIP                } from '../../modules/nf-core/gunzip/main'
include { WINDOWMASKER_MKCOUNTS } from '../../modules/nf-core/windowmasker/mkcounts/main'
include { WINDOWMASKER_USTAT    } from '../../modules/nf-core/windowmasker/ustat/main'


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
    ch_versions = ch_versions.mix ( GUNZIP.out.versions.first() )


    //
    // LOGIC: Extract the genome size for decision making downstream
    //
    ch_genomes_for_gunzip.skip
    | mix( GUNZIP.out.gunzip )
    | map { meta, fa -> [ meta + [genome_size: fa.size()], fa] }
    | set { ch_genome }


    //
    // MODULES: Mask the genome if needed
    //
    if ( params.mask ) {
        WINDOWMASKER_MKCOUNTS ( ch_genome )
        ch_versions = ch_versions.mix ( WINDOWMASKER_MKCOUNTS.out.versions.first() )

        WINDOWMASKER_USTAT ( WINDOWMASKER_MKCOUNTS.out.counts, ch_genome )
        ch_versions = ch_versions.mix ( WINDOWMASKER_USTAT.out.versions.first() )

        ch_fasta = WINDOWMASKER_USTAT.out.intervals
    } else {
        ch_fasta = ch_genome
    }


    emit:
    genome   = ch_fasta            // channel: [ meta, path(genome) ]
    versions = ch_versions         // channel: [ versions.yml ]
}
