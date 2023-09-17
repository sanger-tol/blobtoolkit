// 
// Optional alignment subworkflow using Minimap2
//

include { SAMTOOLS_FASTA                  } from '../../../modules/nf-core/samtools/fasta/main'
include { MINIMAP2_ALIGN as MINIMAP2_HIC  } from '../../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ILMN } from '../../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_CCS  } from '../../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_CLR  } from '../../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ONT  } from '../../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_SORT                   } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX                  } from '../../../modules/nf-core/samtools/index/main'


workflow MINIMAP2_ALIGNMENT {
    take:
    input      // channel: [ val(meta), path(datafile) ]
    fasta      // channel: [ val(meta), path(fasta) ]


    main:
    ch_versions = Channel.empty()


    // Convert reads to FASTA
    SAMTOOLS_FASTA ( input, true )
    ch_versions = ch_versions.mix(SAMTOOLS_FASTA.out.versions.first())


    // Branch input by sequencing type
    SAMTOOLS_FASTA.out.interleaved
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
    fasta
    | map { meta, genome -> genome }
    | set { ch_ref }

    MINIMAP2_HIC ( ch_input.hic, ch_ref, true, false, false )
    ch_versions = ch_versions.mix(MINIMAP2_HIC.out.versions.first())
    
    MINIMAP2_ILMN ( ch_input.illumina, ch_ref, true, false, false )
    ch_versions = ch_versions.mix(MINIMAP2_ILMN.out.versions.first())

    MINIMAP2_CCS ( ch_input.pacbio, ch_ref, true, false, false )
    ch_versions = ch_versions.mix(MINIMAP2_CCS.out.versions.first())

    MINIMAP2_CLR ( ch_input.clr, ch_ref, true, false, false )
    ch_versions = ch_versions.mix(MINIMAP2_CLR.out.versions.first())

    MINIMAP2_ONT ( ch_input.ont, ch_ref, true, false, false )
    ch_versions = ch_versions.mix(MINIMAP2_ONT.out.versions.first())


    // Index aligned reads
    Channel.empty()
    | mix ( MINIMAP2_HIC.out.bam )
    | mix ( MINIMAP2_ILMN.out.bam )
    | mix ( MINIMAP2_CCS.out.bam )
    | mix ( MINIMAP2_CLR.out.bam )
    | mix ( MINIMAP2_ONT.out.bam )
    | set { ch_aligned }

    SAMTOOLS_INDEX ( ch_aligned )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())


    // Combine aligned reads and indices
    ch_aligned
    | join ( SAMTOOLS_INDEX.out.csi )
    | set { bam_csi }


    emit:
    bam_csi                      // channel: [ val(meta), bam, csi ]
    versions = ch_versions       // channel: [ versions.yml ]
}
