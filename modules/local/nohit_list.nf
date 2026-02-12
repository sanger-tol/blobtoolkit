process NOHIT_LIST {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
    'biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(blast)      //path to blast output table in txt format
    tuple val(meta2), path(fasta)      //path to genome fasta file

    output:
    tuple val(meta), path ('*.nohit.txt') , emit: nohitlist
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in sanger-tol/blobtoolkit/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    nohitlist.sh ${fasta} ${blast} ${prefix} $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nohit_list: 1.0
    END_VERSIONS
    """
}
