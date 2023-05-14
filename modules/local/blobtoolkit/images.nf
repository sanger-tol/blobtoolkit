process BLOBTOOLKIT_IMAGES {
    tag "${meta.id}_${plot}"
    label 'process_single'

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "BLOBTOOLKIT_IMAGES module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    container "genomehubs/blobtoolkit:4.1.5"

    input:
    tuple val(meta), path(blobdir)
    each plot

    output:
    tuple val(meta), path('*.png') , emit: png
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    blobtk plot \\
        -v ${plot} \\
        -o ${prefix}.${plot}.png \\
        -d ${blobdir} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: \$(btk --version | cut -d' ' -f2 | sed 's/v//')
    END_VERSIONS
    """
}
