process NOHIT_LIST {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
    'quay.io/biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(blast)      //path to blast output table in txt format
    tuple val(meta), path(fasta)      //path to genome fasta file

    output:
    tuple val(meta), path ('*.nohit.txt') , emit: nohitlist
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    grep '>' ${fasta} | \\
        grep -v -w -f <(awk -v evalue="$args" '{{if($14<{evalue}){{print $1}}}}' ${blast} | sort | uniq) | \\
        cut -f1 | sed 's/>//' > ${prefix}.nohit.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nohit_list: 1.0
    END_VERSIONS
    """
}
