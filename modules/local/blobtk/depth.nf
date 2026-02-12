process BLOBTK_DEPTH {
    tag "${meta.id}"
    label 'process_single'

    container "docker.io/genomehubs/blobtk:0.5.1"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path('*.regions.bed.gz') , emit: bed
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "BLOBTOOLKIT_DEPTH module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    blobtk depth \\
        -b ${bam} \\
        $args \\
        -O ${prefix}.regions.bed.gz \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtk: \$(blobtk --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
