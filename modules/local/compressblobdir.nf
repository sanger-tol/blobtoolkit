process COMPRESSBLOBDIR {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::pigz=2.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.8':
        'biocontainers/pigz:2.8' }"

    input:
    tuple val(meta), path(input, stageAs: "input_blobdir")
    tuple val(meta1), path(summary_json)
    tuple val(meta2), path(meta_json)

    output:
    tuple val(meta), path(prefix), emit: blobdir
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    cp ${input}/* ${prefix}/
    cp ${summary_json} ${prefix}/summary.json
    cp ${meta_json} ${prefix}/meta.json
    pigz --processes $task.cpus ${prefix}/*.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz:\$(echo \$(pigz --version 2>&1) | sed 's/^.*pigz\\w*//' )
    END_VERSIONS
    """
}
