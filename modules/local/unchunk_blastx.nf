process UNCHUNK_BLASTX {
    tag "$fasta"

    container "genomehubs/blobtoolkit-blobtools"

    input:
    path raw_proteomes

    output:
    path "reference_proteomes.out" , emit: proteomes
    path "versions.yml"            , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    btk pipeline unchunk-blast \\
        $args \\
        --in ${raw_proteomes} \\
        --out reference_proteomes.out 2> unchunk_blastx.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: \$( echo "3.1.0" | sed 's/blobtoolkit //g')
    END_VERSIONS
    """
}
