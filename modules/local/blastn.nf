process BLASTN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::blast=2.15.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.15.0--pl5321h6f7f691_1' :
        'biocontainers/blast:2.15.0--pl5321h6f7f691_1' }"

    input:
    tuple val(meta), path(fasta)
    path  db
    val   taxid

    output:
    tuple val(meta), file('*.blastn.txt'), emit: txt
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def exclude_taxon = taxid ? "-negative_taxids ${taxid}" : ''
    """
    DB=`find -L ./ -name "*.ndb" | sed 's/\\.ndb\$//'`
    blastn \\
        -num_threads $task.cpus \\
        -db \$DB \\
        -query $fasta \\
        $exclude_taxon \\
        $args \\
        -out ${prefix}.blastn.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
