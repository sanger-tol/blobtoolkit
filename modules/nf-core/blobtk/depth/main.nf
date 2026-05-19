process BLOBTK_DEPTH {
    tag "${meta.id}"
    label 'process_single'

    container "docker.io/genomehubs/blobtk:0.8.1"

    input:
    tuple val(meta), path(bam), path(index)

    output:
    tuple val(meta), path('*.regions.bed.gz') , emit: bed
    tuple val("${task.process}"), val("blobtk"), eval("blobtk --version | cut -d' ' -f2"), topic: versions, emit: versions_blobtk

    when:
    task.ext.when == null || task.ext.when

    script:

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "BLOBTK_DEPTH module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    """
    blobtk depth \\
        -b ${bam} \\
        $args \\
        -O ${prefix}.regions.bed.gz \\
    """

    stub:
    def prefix      = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.regions.bed.gz
    """
}
