process UNCHUNK_BLASTX {
    tag "$fasta"

    container "genomehubs/blobtoolkit-blobtools"

    input:
    path raw_blastn
    val max_target_seqs

    output:
    path "*.blastn"     , emit: proteomes
    path "versions.yml" , emit: versions

    script:
    def fname  = raw_blastn.simpleName()
    def log    = "${fname}.log"
    def output = "${fname}.blastn"
    """
    btk pipeline unchunk-blast \\
        --in ${raw_blastn} \\
        --count $max_target_seqs \\
        --out ${output} 2> ${log}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: \$( echo "3.1.0" | sed 's/blobtoolkit //g')
    END_VERSIONS
    """
}
