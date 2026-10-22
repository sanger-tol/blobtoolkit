process BLOBTOOLKIT_SUMMARY {
    tag "${meta.id}"
    label 'process_single'

    container "docker.io/genomehubs/blobtoolkit:4.4.6"

    input:
    tuple val(meta), path(blobdir)

    output:
    tuple val(meta), path("*.json"), emit: json
    tuple val("${task.process}"), val("blobtoolkit"), eval("btk --version | cut -d' ' -f2 | sed 's/v//'"), topic: versions, emit: versions_blobtoolkit

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "BLOBTOOLKIT_SUMMARY module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    blobtools filter \\
        ${args} \\
        --summary ${prefix}.summary.json ${blobdir}
    """
}
