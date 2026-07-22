process CREATE_BED {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
    'biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(tsv)      //path to tsv output from fasta windows

    output:
    tuple val(meta), path ('*.bed') , emit: bed
    tuple val("${task.process}"), val("createbed"), val("1.03"), topic: versions, emit: versions_createbed

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cut -f 1,2,3 $tsv | sed '1d' > ${prefix}.bed
    """
}
