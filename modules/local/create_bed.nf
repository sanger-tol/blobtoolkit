process CREATE_BED {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "conda-forge::gawk=5.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
    'quay.io/biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(tsv)      //path to tsv output from fasta windows

    output:
    path '*.bed'       , emit: bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cut -f 1,2,3 $tsv | sed '1d' $args > ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        create_bed: 1.01
    END_VERSIONS
    """
}