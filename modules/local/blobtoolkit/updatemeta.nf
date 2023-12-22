process BLOBTOOLKIT_UPDATEMETA {
    tag "$meta.id"
    label 'process_single'

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "BLOBTOOLKIT_UPDATEMETA module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    container "docker.io/pacificbiosciences/pyyaml:5.3.1"

    input:
    tuple val(meta), path(input)
    path versions

    output:
    tuple val(meta), path(prefix), emit: blobdir
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    update_versions.py \\
        ${args} \\
        --meta ${input}/meta.json \\
        --software ${versions} \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        update_versions.py: \$(update_versions.py --version | cut -d' ' -f2)
    END_VERSIONS
    """

}
