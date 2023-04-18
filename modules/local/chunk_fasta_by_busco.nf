process CHUNK_FASTA_BY_BUSCO {
    tag "$meta.id"
    label 'process_single'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the CHUNK_FASTA_BY_BUSCO module. Please use docker or singularity containers."
    }
    container 'genomehubs/blobtoolkit:4.1.2'

    input:
    tuple val(meta), path(fasta)
    path busco_table

    output:
    tuple val(meta), path('*.output_chunks.fasta') , emit: chunks
    path "versions.yml"                            , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    btk pipeline chunk-fasta \\
        $args \\
        --in ${fasta} \\
        --busco ${busco_table} \\
        --out ${prefix}.chunks.fasta \\
        --bed None
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
         blobtoolkit: \$(btk --version | cut -d' ' -f2 | sed 's/v//')
    END_VERSIONS
    """
}
