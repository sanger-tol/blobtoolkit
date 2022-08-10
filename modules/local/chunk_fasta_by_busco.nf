process CHUNK_FASTA_BUSCO {
    tag "$fasta"

    container "genomehubs/blobtoolkit-blobtools"

    input:
    path fasta
    path busco_table
    val chunk
    val overlap
    val max_chunks
    val min_length

    output:
    path "output_chunks.fasta" , emit: chunks
    path "versions.yml"        , emit: versions

    script:
    """
    btk pipeline chunk-fasta \\
        --in ${fasta} \\
        --busco ${busco_table} \\
        --chunk $chunk \\
        --overlap $overlap \\
        --max-chunks $max_chunks \\
        --min-length $min_length \\
        --out output.chunks.fasta \\
        --bed None 2> chunk_fasta.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: \$( echo "3.1.0" | sed 's/blobtoolkit //g')
    END_VERSIONS
    """
}
