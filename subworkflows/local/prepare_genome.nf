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
        ch_unzipped = GUNZIP ( fasta ).gunzip
        ch_versions = ch_versions.mix ( GUNZIP.out.versions )
    } else {
        ch_unzipped = fasta
    }

    //
    // LOGIC: Extract the genome size for decision making downstream
    //
    ch_unzipped
    | map { meta, fa -> [ meta + [genome_size: fa.size()], fa] }
    | set { ch_genome }

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
