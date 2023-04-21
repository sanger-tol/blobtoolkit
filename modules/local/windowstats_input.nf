process WINDOWSTATS_INPUT {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2':
        'quay.io/biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(freq)
    tuple val(meta), path(mononuc)
    tuple val(meta), path(mosdepth)
    tuple val(meta), path(countbusco)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    windowstats_input.py \\
        --freq ${freq} \\
        --mononuc ${mononuc} \\
        --mosdepth ${mosdepth} \\
        --countbusco ${countbusco} \\
        --output ${prefix}.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        windowstats_input.py: \$(windowstats_input.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
