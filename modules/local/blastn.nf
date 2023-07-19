process BLASTN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::blast=2.13.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.13.0--hf3cf87c_0' :
        'biocontainers/blast:2.13.0--hf3cf87c_0' }"

    input:
    tuple val(meta), path(fasta)
    path  db
    val   taxid

    output:
    tuple val(meta), path('*.blastn.txt'), emit: txt
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = "${prefix}.blastn.txt"
    def log = "${prefix}.blastn.log"
    """
    DB=`find -L ./ -name "*.ndb" | sed 's/\\.ndb\$//'`
    if [ -s $fasta ]; then \\
        blastn \\
            -num_threads $task.cpus \\
            -db \$DB \\
            -query $fasta \\
            -negative_taxids $taxid \\
            $args \\
            > $output 2> $log || \\
            sleep 30; \\
            if [ -s $log ]; then \\
                echo "Restarting blastn without taxid filter" >> $log; \\
                > $output; \\
                blastn \\
                    -num_threads $task.cpus \\
                    -db \$DB \\
                    -query $fasta \\
                    $args \\
                    > $output 2>> $log; \\
            fi \\
            else \\
            > $output; \\
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
