process EXTRACT_BUSCO_GENES {
    tag "$meta.id"

    container "genomehubs/blobtoolkit-blobtools:3.3.4"

    input:
    tuple val(meta), val(tables)

    output:
    tuple val(meta), path('*.busco_genes.fasta') , emit: fasta
    path "versions.yml"                          , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    btk pipeline extract-busco-genes \\
        --busco $tables \\
        --out ${prefix}_busco_genes.fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: \$( echo "3.1.0" | sed 's/blobtoolkit //g')
    END_VERSIONS
    """
}