process TARGZ {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::tar=1.34"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.tar.gz"), emit: archive
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${input}"
    """
    tar -czvf ${prefix}.tar.gz $input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tar: \$(echo \$(tar --version 2>&1) | sed 's/^.*(GNU tar) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
