process EXTRACT_BUSCO_GENES {
    tag "$fasta"

    container "genomehubs/blobtoolkit-blobtools"

    input:
    path dir

    output:
    path "output_busco_genes.fasta" , emit: fasta
    path "versions.yml"             , emit: versions

    script:
    def tables = Channel.fromPath( ["$dir/**/run_archaea_odb10/full_table.tsv", "$dir/**/run_bacteria_odb10/full_table.tsv", "$dir/**/run_eukaryota_odb10/full_table.tsv"] )
    def busco_tables = tables.toList()
    """
    btk pipeline extract-busco-genes \\
        --busco $busco_tables \\
        --out output_busco_genes.fasta 2> extract_busco_genes.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: \$( echo "3.1.0" | sed 's/blobtoolkit //g')
    END_VERSIONS
    """
}
