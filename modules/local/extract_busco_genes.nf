process EXTRACT_BUSCO_GENES {
    tag "$fasta"

    container "genomehubs/blobtoolkit-blobtools"

    input:
    path busco_table

    output:
    path 'output_busco_genes.fasta' , emit: fasta
    path "versions.yml"             , emit: versions

    script:
    // busco_table input is aggregated using a Snakemake function called expand():
    // see input in https://github.com/blobtoolkit/blobtoolkit/blob/main/src/btk_pipeline/rules/extract_busco_genes.smk
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
