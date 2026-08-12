//
// Generate summary and static plots from blobdir
//

include { BLOBTOOLKIT_SUMMARY } from '../../modules/local/blobtoolkit/summary'
include { BLOBTK_PLOT         } from '../../modules/nf-core/blobtk/plot/main'

workflow VIEW {
    take:
    blobdir     // channel: [ val(meta), path(blobdir) ]


    main:
    //
    // MODULE: GENERATE SUMMARY FILE
    //
    BLOBTOOLKIT_SUMMARY ( blobdir )


    //
    // MODULE: GENERATE STATIC PLOTS IN PNG/SVG FORMAT
    //
    plots = channel.of(
        [
            name: "blob",
            args: "-v blob"
        ],
        [
            name: "grid",
            args: "-v blob --shape grid -x position"
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
        // Only add "-w 0.01" to the grid paramters if there are such windows available
        .map { meta, local, btk_args ->
            ((btk_args.name == "grid") && file(local.toUriString() + "/length_windows_0.01.json").exists())
            ? [meta, local, [name: "grid", args: btk_args.args + " -w 0.01"]]
            : [meta, local, btk_args]
        }
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
}
