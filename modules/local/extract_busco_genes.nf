process EXTRACT_BUSCO_GENES {
    tag "$busco_table"

    container "${ workflow.containerEngine == 'docker' ?
        'genomehubs/blobtoolkit-blobtools' :
        '' }"

    input:
    path busco_table

    output:
    path 'output_busco_genes.fasta' , emit: fasta
    path "versions.yml"             , emit: versions

    script:
    // busco_table might require formatting
    // This script is bundled with the pipeline, in nf-core/blobtoolkit/bin/
    """
    btk pipeline extract-busco-genes \\
        --busco $busco_table \\
        --out output_busco_genes.fasta 2> extract_busco_genes.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: \$( echo "3.1.0" | sed 's/blobtoolkit //g')
    END_VERSIONS
    """
}
