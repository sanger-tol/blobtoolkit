process GET_NOHIT_LIST {
    tag "$meta.id"

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta), path(blastx)
    val  evalue

    output:
    path "*.nohit.txt"  , emit: nohit_list
    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    grep '>' ${fasta} | \\
        grep -v -w -f <(awk '{{if($14<{$evalue}){{print $1}}}}' ${blastx} | \\
        sort | uniq) | \\
        cut -f1 | sed 's/>//' > ${prefix}.nohit.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ubuntu: \$( echo "20.04" | sed 's/ubuntu //g')
    END_VERSIONS
    """
}
