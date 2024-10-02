//
// Optional alignment subworkflow using Minimap2
//

include { MINIMAP2_ALIGN as MINIMAP2_HIC  } from '../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ILMN } from '../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_CCS  } from '../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_CLR  } from '../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ONT  } from '../../modules/nf-core/minimap2/align/main'


workflow MINIMAP2_ALIGNMENT {
    take:
    input      // channel: [ val(meta), path(datafile) ]
    fasta      // channel: [ val(meta), path(fasta) ]


    main:
    ch_versions = Channel.empty()


    // Branch input by sequencing type
    input
    | branch {
        meta, reads ->
            hic: meta.datatype == "hic"
            illumina : meta.datatype == "illumina"
            pacbio : meta.datatype == "pacbio"
            clr : meta.datatype == "pacbio_clr"
            ont : meta.datatype == "ont"
    }
    | set { ch_input }


    // Align with Minimap2
    MINIMAP2_HIC ( ch_input.hic, fasta, true, false, false )
    ch_versions = ch_versions.mix(MINIMAP2_HIC.out.versions.first())

    MINIMAP2_ILMN ( ch_input.illumina, fasta, true, false, false )
    ch_versions = ch_versions.mix(MINIMAP2_ILMN.out.versions.first())

    MINIMAP2_CCS ( ch_input.pacbio, fasta, true, false, false )
    ch_versions = ch_versions.mix(MINIMAP2_CCS.out.versions.first())

    MINIMAP2_CLR ( ch_input.clr, fasta, true, false, false )
    ch_versions = ch_versions.mix(MINIMAP2_CLR.out.versions.first())

    MINIMAP2_ONT ( ch_input.ont, fasta, true, false, false )
    ch_versions = ch_versions.mix(MINIMAP2_ONT.out.versions.first())


    // Combine aligned reads
    Channel.empty()
    | mix ( MINIMAP2_HIC.out.bam )
    | mix ( MINIMAP2_ILMN.out.bam )
    | mix ( MINIMAP2_CCS.out.bam )
    | mix ( MINIMAP2_CLR.out.bam )
    | mix ( MINIMAP2_ONT.out.bam )
    | set { ch_aligned }


    emit:
    aln      = ch_aligned        // channel: [ val(meta), bam ]
    versions = ch_versions       // channel: [ versions.yml ]
}
