process BLOBTK_IMAGES {
    tag "${meta.id}_${plot}"
    label 'process_single'

    container "docker.io/genomehubs/blobtk:0.5.1"

    input:
    tuple val(meta), path(blobdir)
    each plot
    val format

    output:
    tuple val(meta), path('*.png') , optional: true, emit: png
    tuple val(meta), path('*.svg') , optional: true, emit: svg
    tuple val("${task.process}"), val("blobtk"), eval("blobtk --version | cut -d' ' -f2"), topic: versions, emit: versions_blobtk

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "BLOBTOOLKIT_IMAGES module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def legend = plot.equals("snail") ? "" : "--legend full"
    """
    blobtk plot \\
        -v ${plot} \\
        -d ${blobdir} \\
        -o ${prefix}.${plot}.${format} \\
        ${legend} \\
        $args
    """
}
