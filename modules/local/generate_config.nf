process GENERATE_CONFIG {
    tag "$meta.id"
    label 'process_single'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the GENERATE_CONFIG module. Please use docker or singularity containers."
    }
    container 'genomehubs/blobtoolkit:4.0.7'

    input:
    val meta

    output:
    tuple val(meta), path("${prefix}/*.yaml") , emit: yaml
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    blobtoolkit-pipeline generate-config ${prefix}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit-pipeline: \$(blobtoolkit-pipeline --version | cut -d' ' -f2 | sed 's/v//')
    END_VERSIONS
    """
}
