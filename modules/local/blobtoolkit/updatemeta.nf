process BLOBTOOLKIT_UPDATEMETA {
    tag "$meta.id"
    label 'process_single'

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "BLOBTOOLKIT_UPDATEMETA module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    container "docker.io/genomehubs/blobtoolkit:4.3.9"

    input:
    tuple val(meta), path(input)
    path versions

    output:
    tuple val(meta), path("*.json"), emit: json
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    update_versions.py \\
        ${args} \\
        --meta_in ${input}/meta.json \\
        --software ${versions} \\
        --meta_out ${prefix}.meta.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        update_versions.py: \$(update_versions.py --version | cut -d' ' -f2)
    END_VERSIONS
    """

}
