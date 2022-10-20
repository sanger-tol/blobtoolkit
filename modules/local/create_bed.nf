process CREATE_BED {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "conda-forge::python=3.9.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path (fasta)      //path to fasta windows file

    output:
    path '*.bed'       , emit: bed
    //path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cut -f 1,2,3 $fasta | sed '1d' >> ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        create_bed 1.0
    END_VERSIONS
    """
}