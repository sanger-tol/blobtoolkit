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

    output:
    tuple val(meta), path("*.yaml"), emit: yaml
    tuple val(meta), path("*.csv") , emit: csv
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def busco_param = busco_lin ? "--busco '${busco_lin}'" : ""
    """
    generate_config.py $fasta "$taxon_query" $lineage_tax_ids $busco_param ${prefix}.yaml ${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        generate_config.py: \$(generate_config.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
