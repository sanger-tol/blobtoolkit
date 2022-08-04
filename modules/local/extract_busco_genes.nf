process EXTRACT_BUSCO_GENES {
    tag "$fasta"

    container "genomehubs/blobtoolkit-blobtools"

    input:
    val tables

    output:
    path "output_busco_genes.fasta" , emit: fasta
    path "versions.yml"             , emit: versions

    script:
    """
    btk pipeline extract-busco-genes \\
        --busco $tables \\
        --out output_busco_genes.fasta 2> extract_busco_genes.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: \$( echo "3.1.0" | sed 's/blobtoolkit //g')
    END_VERSIONS
    """
}
