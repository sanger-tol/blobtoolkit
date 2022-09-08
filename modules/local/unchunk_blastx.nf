process UNCHUNK_BLASTX {
    tag "$fasta"

    container "genomehubs/blobtoolkit-blobtools:3.3.4"

    input:
    val prefix
    path raw_proteomes

    output:
    path "*.reference_proteomes.out" , emit: proteomes
    path "versions.yml"              , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    btk pipeline unchunk-blast \\
        $args \\
        --in ${raw_proteomes} \\
        --out ${prefix}.reference_proteomes.out 2> ${prefix}.unchunk_blastx.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
         blobtoolkit: \$(btk --version | cut -d' ' -f2 | sed 's/v//')
    END_VERSIONS
    """
}
