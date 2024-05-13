process NCBI_GET_ODB_TAXON {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::requests=2.26.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/requests:2.26.0':
        'biocontainers/requests:2.26.0' }"

    input:
    tuple val(meta), val(taxon_query)
    val busco_lin
    path lineage_tax_ids

    output:
    tuple val(meta), path("*.csv"), emit: csv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def busco_param = busco_lin ? "--busco '${busco_lin}'" : ""
    """
    get_odb_taxon.py "$taxon_query" $lineage_tax_ids $busco_param ${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        get_odb_taxon.py: \$(get_odb.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
