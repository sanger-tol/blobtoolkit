    process RUN_BLASTN {
    tag "$fasta"

    container "quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0"

    input:
    path fasta
    path db
    val evalue
    val max_target_seqs

    output:
    path "*.blastn.raw"  , emit: blastn_out
    path "versions.yml" , emit: versions

    script:
    def fname  = fasta.simpleName()
    def log    = "${fname}.log"
    def output = "${fname}.blastn.raw"
    """
    if [ -s ${fasta} ]; then \\
            blastn -task megablast \\
                -query ${fasta} \\
                -db ${db} \\
                -outfmt "6 qseqid staxids bitscore std" \\
                -max_target_seqs $max_target_seqs \\
                -max_hsps 1 \\
                -evalue $evalue \\
                -num_threads $task.cpus -negative_taxids \\
                -lcase_masking \\
                -dust "20 64 1" \\
                > ${output} 2> ${log} || \\
            sleep 30; \
            if [ -s ${log} ]; then \\
                echo "Restarting blastn without taxid filter" >> ${log}; \\
                > ${output}; \
                blastn -task megablast \\
                    -query ${fasta} \\
                    -db ${db} \\
                    -outfmt "6 qseqid staxids bitscore std" \\
                    -max_target_seqs $max_target_seqs \\
                    -max_hsps 1 \\
                    -evalue $evalue \\
                    -num_threads $task.cpus \\
                    -lcase_masking \\
                    -dust "20 64 1" \\
                    > ${output} 2>> ${log}; \\
            fi \\
        else \\
            > ${output}; \\
        fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
