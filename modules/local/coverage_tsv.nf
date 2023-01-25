process COVERAGE_TSV {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    tuple val(meta), path(mosdepth)         //taking output of mosdepth
    tuple val(meta), path(countbusogenes)   //taking output of count_buscogenes

    output:
    tuple val(meta), path ('*_coverage.tsv')       , emit: cov_tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    sed -i 1i"Sequence  Start   End ${prefix}_cov" ${mosdepth}
    paste ${countbusogenes} ${mosdepth} > ${prefix}_coverage.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coverage_tsv : 1.00
    END_VERSIONS
    """
}