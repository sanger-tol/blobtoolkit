process CHUNK_FASTA_BUSCO {
    tag "$fasta"

    container "genomehubs/blobtoolkit-blobtools"

    input:
    path fasta
    path busco_table

    output:
    path "output_chunks.fasta" , emit: chunks
    path "versions.yml"        , emit: versions

    script:
    """
    btk pipeline chunk-fasta \\
        --in ${fasta} \\
        --chunk "${params.chunk}"\\
        --overlap "${params.overlap}" \\
        --max-chunks "${params.max_chunks}" \\
        --min-length "${params.min_length}" \\
        --busco ${busco_table} \\
        --out output.chunks.fasta \\
        --bed None 2> chunk_fasta.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: \$( echo "3.1.0" | sed 's/blobtoolkit //g')
    END_VERSIONS
    """
}
