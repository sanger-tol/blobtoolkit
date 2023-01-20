process EXTRACT_BUSCO_GENES {
    tag "$meta.id"
    label 'process_single'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the EXTRACT_BUSCO_GENES module. Please use docker or singularity containers."
    }
    container "genomehubs/blobtoolkit-blobtools:3.4.2"

    input:
    tuple val(meta), path(arc), path(bac), path(euk)

    output:
    tuple val(meta), path('*_busco_genes.fasta') , emit: fasta
    path "versions.yml"                          , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    btk pipeline extract-busco-genes \\
        --busco $arc/full_table.tsv.gz \\
        --busco $bac/full_table.tsv.gz \\
        --busco $euk/full_table.tsv.gz \\
        --out ${prefix}_busco_genes.fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: \$(btk --version | cut -d' ' -f2 | sed 's/v//')
    END_VERSIONS
    """
}
