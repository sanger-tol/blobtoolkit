process ADD_SUMMARY_TO_METADATA {
    tag "$meta.id"
    label 'process_single'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the ADD_SUMMARY_TO_METADATA module. Please use docker or singularity containers."
    }
    container 'genomehubs/blobtoolkit:4.0.7'

    input:
    tuple val(meta), path(yaml)

    output:
    tuple val(meta), path('*.config.yaml*') , emit: yaml
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    btk pipeline add-summary-to-metadata \\
        --config ${yaml} \\
        --out ${prefix}.config.yaml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit-pipeline: \$(blobtoolkit-pipeline --version | cut -d' ' -f2 | sed 's/v//')
    END_VERSIONS
    """
}
