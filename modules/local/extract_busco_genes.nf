process EXTRACT_BUSCO_GENES {
    tag "$meta.id"

    container "genomehubs/blobtoolkit-blobtools:3.3.4"

    input:
    tuple val(meta), path(arc), path(bac), path(euk)

    output:
    tuple val(meta), path('*.busco_genes.fasta') , emit: fasta
    path "versions.yml"                          , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def arct = "arc_full_table.tsv"
    def bact = "bac_full_table.tsv"
    def eukt = "euk_full_table.tsv"
    def tables = ["\"$arct\"", "\"$bact\"", "\"$eukt\""]
    """
    cp ${arc} $arct
    cp ${bac} $bact
    cp ${euk} $eukt

    btk pipeline extract-busco-genes \\
        --busco $tables \\
        --out ${prefix}_busco_genes.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: \$(btk --version | cut -d' ' -f2 | sed 's/v//')
    END_VERSIONS
    """
}
