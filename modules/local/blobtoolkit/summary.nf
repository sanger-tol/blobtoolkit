process BLOBTOOLKIT_SUMMARY {
    tag "$meta.id"
    label 'process_single'

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "BLOBTOOLKIT_SUMMARY module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    container "docker.io/genomehubs/blobtoolkit:4.4.6"

    input:
    tuple val(meta), path(blobdir)

    output:
    tuple val(meta), path("*.json"), emit: json
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    blobtools filter \\
        ${args} \\
        --summary ${prefix}.summary.json ${blobdir}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: \$(btk --version | cut -d' ' -f2 | sed 's/v//')
    END_VERSIONS
    """
}
