//
// Prepare genome for downstream processing
//

include { GUNZIP                } from '../../modules/nf-core/gunzip/main'
include { WINDOWMASKER_MKCOUNTS } from '../../modules/nf-core/windowmasker/mkcounts/main'
include { WINDOWMASKER_USTAT    } from '../../modules/nf-core/windowmasker/ustat/main'


workflow PREPARE_GENOME {
    take:
    fasta     // channel: [ meta, path(genome) ]


    main:
    ch_versions = Channel.empty()


    //
    // MODULE: Decompress FASTA file if needed
    //
    if ( params.fasta.endsWith('.gz') ) {
        ch_genome   = GUNZIP ( fasta ).gunzip
        ch_versions = ch_versions.mix ( GUNZIP.out.versions )
    } else {
        ch_genome   = fasta
    }


    //
    // MODULES: Mask the genome if needed
    //
    if ( params.mask ) {
        WINDOWMASKER_MKCOUNTS ( ch_genome )
        ch_versions = ch_versions.mix ( WINDOWMASKER_MKCOUNTS.out.versions )

        WINDOWMASKER_USTAT ( WINDOWMASKER_MKCOUNTS.out.counts, ch_genome )
        ch_versions = ch_versions.mix ( WINDOWMASKER_USTAT.out.versions )

        ch_fasta = WINDOWMASKER_USTAT.out.intervals
    } else {
        ch_fasta = ch_genome
    }

    
    emit:
    genome   = ch_fasta            // channel: [ meta, path(genome) ]
    versions = ch_versions         // channel: [ versions.yml ]
}