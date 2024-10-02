process GENERATE_CONFIG {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::requests=2.28.1 conda-forge::pyyaml=6.0"
    container "docker.io/genomehubs/blobtoolkit:4.3.9"

    input:
    tuple val(meta), val(fasta)
    val taxon_query
    val busco_lin
    path lineage_tax_ids
    tuple val(meta2), path(blastn)
    val reads
    // The following are passed as "val" because we just want to know the full paths. No staging necessary
    val blastp_path
    val blastx_path
    val blastn_path
    val taxdump_path

    output:
    tuple val(meta), path("*.yaml")          , emit: yaml
    tuple val(meta), path("*.csv")           , emit: csv
    tuple val(meta), path("*.synonyms.tsv")  , emit: synonyms_tsv,   optional: true
    tuple val(meta), path("*.categories.tsv"), emit: categories_tsv, optional: true
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def busco_param = busco_lin ? "--busco '${busco_lin}'" : ""
    def accession_params = params.accession ? "--accession ${params.accession}" : ""
    def input_reads = reads.collect{"--read_id ${it[0].id} --read_type ${it[0].datatype} --read_layout ${it[0].layout} --read_path ${it[1]}"}.join(' ')
    """
    generate_config.py \\
        --fasta $fasta \\
        --taxon_query "$taxon_query" \\
        --lineage_tax_ids $lineage_tax_ids \\
        $busco_param \\
        $accession_params \\
        --nt $blastn \\
        $input_reads \\
        --blastp ${blastp_path} \\
        --blastx ${blastx_path} \\
        --blastn ${blastn_path} \\
        --taxdump ${taxdump_path} \\
        --output_prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        generate_config.py: \$(generate_config.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
