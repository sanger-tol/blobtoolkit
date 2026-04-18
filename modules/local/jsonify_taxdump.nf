process JSONIFY_TAXDUMP {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::requests=2.28.1 conda-forge::pyyaml=6.0"
    container "docker.io/genomehubs/blobtoolkit:4.4.6"

    input:
    tuple val(meta), path(taxdump)

    output:
    tuple val(meta), path("*.json") , emit: json
    tuple val("${task.process}"), val('jsonify_taxdump.py'), eval('jsonify_taxdump.py --version 2>&1 | cut -d" " -f2'), emit: versions_jsonify_taxdump, topic: versions

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
    """
}
