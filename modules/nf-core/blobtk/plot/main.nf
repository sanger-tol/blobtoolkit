process BLOBTK_PLOT {
    tag "$prefix"
    label 'process_single'

    container "docker.io/genomehubs/blobtk:0.8.1"

    input:
    tuple val(meta), path(fasta)
    path(local_path)                // Genuine path location must be a path.
    val(online_path)                // HTTPS location needs to remain a value
    val extra_args                  // In format [name: "", args: ""]
    val format                      // Output format, e.g. png or svg

    output:
    tuple val(meta), path("*.png"), optional: true, emit: png
    tuple val(meta), path("*.svg"), optional: true, emit: svg
    tuple val("${task.process}"), val("blobtk"), eval("blobtk --version | cut -d' ' -f2"), topic: versions, emit: versions_blobtk

    when:
    task.ext.when == null || task.ext.when

    script:

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "BLOBTK_PLOT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args    = task.ext.args ?: ''

    if ( online_path && local_path ) {
        error "BLOBTK_PLOT can't use both local_path and online_path, use `[]` as input for the unused channel."
    }

    def resource = online_path ?: local_path
    def legend   = extra_args.args.contains("-v snail") ? "" : "--legend full"

    prefix       = task.ext.prefix ?: "${meta.id}"

    """
    blobtk plot \\
        -d ${resource} \\
        -o ${prefix}.${format} \\
        ${legend} \\
        $args
    """

    stub:
    prefix      = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.${format}
    """
}
