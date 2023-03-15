process GENERATE_SUMMARY {
    tag "$meta.id"
    label 'process_single'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the GENERATE_SUMMARY module. Please use docker or singularity containers."
    }
    container 'genomehubs/blobtoolkit:4.0.7'

    input:
    tuple val(meta), path(blobdir)
    tuple val(meta), path(png)

    output:
    tuple val(meta), path("${blobdir}/summary.json") , emit: json
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def output = "${blobdir}/summary.json"
    """
    blobtools filter \\
        --summary $output ${blobdir} \\
        $args
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit-pipeline: \$(blobtoolkit-pipeline --version | cut -d' ' -f2 | sed 's/v//')
    END_VERSIONS
    """
}