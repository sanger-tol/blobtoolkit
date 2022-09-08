process CHUNK_FASTA_BUSCO {
    tag "$fasta"

    container "genomehubs/blobtoolkit-blobtools:3.3.4"

    input:
    val  prefix
    path fasta
    path busco_table

    output:
    path "*.output_chunks.fasta" , emit: chunks
    path "versions.yml"          , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    btk pipeline chunk-fasta \\
        $args \\
        --in ${fasta} \\
        --busco ${busco_table} \\
        --out ${prefix}.chunks.fasta \\
        --bed None 2> ${prefix}.chunk_fasta.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
         blobtoolkit: \$(btk --version | cut -d' ' -f2 | sed 's/v//')
    END_VERSIONS
    """
}
