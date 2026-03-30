//
// Generate summary and static plots from blobdir
//

include { BLOBTOOLKIT_SUMMARY } from '../../modules/local/blobtoolkit/summary'
include { BLOBTK_PLOT         } from '../../modules/nf-core/blobtk/plot/main'

workflow VIEW {
    take:
    blobdir     // channel: [ val(meta), path(blobdir) ]


    main:
    ch_versions = channel.empty()


    //
    // Generate summary file
    //
    BLOBTOOLKIT_SUMMARY ( blobdir )
    ch_versions = ch_versions.mix ( BLOBTOOLKIT_SUMMARY.out.versions.first() )


    //
    // MODULE: Generate static plots in png/svg format
    //
    plots = channel.of(
        [
            name: "blob",
            args: "-v blob"
        ],
        [
            name: "cumulative",
            args: "-v cumulative"
        ],
        [
            name: "snail",
            args: "-v snail"
        ]
    )

    ch_blobtk_plot_input = blobdir
        .combine(plots)
        .multiMap { meta, local, btk_args ->
            fasta: [meta, []]
            local_path: local
            online_path: []
            args: btk_args
        }

    BLOBTK_PLOT(
        ch_blobtk_plot_input.fasta,
        ch_blobtk_plot_input.local_path,
        ch_blobtk_plot_input.online_path,
        ch_blobtk_plot_input.args,
        params.image_format
    )

    ch_images = BLOBTK_PLOT.out.png.mix(BLOBTK_PLOT.out.svg)

    emit:
    summary  = BLOBTOOLKIT_SUMMARY.out.json  // channel: [ val(meta), path(json) ]
    images   = ch_images                     // channel: [ val(meta), path(png/svg) ]
    versions = ch_versions                   // channel: [ versions.yml ]
}
