process JSONIFY_TAXDUMP {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::requests=2.28.1 conda-forge::pyyaml=6.0"
    container "docker.io/genomehubs/blobtoolkit:4.4.4"

    input:
    tuple val(meta), path(taxdump)

    output:
    tuple val(meta), path("*.json") , emit: json
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    jsonify_taxdump.py \\
        $taxdump \\
        $args \\
        > ${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jsonify_taxdump.py: \$(jsonify_taxdump.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
