// Based on the EXTRACT_BUSCO_GENES module 

process COUNT_BUSCOGENES {
    tag "$meta.id"
    label 'process_single'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the COUNT_BUSCOGENES module. Please use docker or singularity containers."
    }
    container "genomehubs/blobtoolkit-blobtools:3.4.2"

    input:
    tuple val(meta), path(tsv) //full table tsv files
    tuple val(meta), path(bed) //modified fasta windows bed file

    output:
    tuple val(meta), path('*_busco_genes_count.tsv') , emit: tsv //output tsv
    path "versions.yml"                              , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def busco_command = tsv.collect{"--in $it"}.join(' ')

    """
    btk pipeline count-busco-genes \\
            $busco_command \\
            --mask ${bed} \\
            --out ${prefix}__busco_genes_count.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        count_buscogenes: 1.00
    END_VERSIONS
    """
}