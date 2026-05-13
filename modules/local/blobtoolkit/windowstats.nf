process BLOBTOOLKIT_WINDOWSTATS {
    tag "${meta.id}"

    container "docker.io/genomehubs/blobtoolkit:4.4.6"

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path('*_window_stats*.tsv') , emit: tsv
    tuple val("${task.process}"), val("blobtoolkit"), eval("btk --version | cut -d' ' -f2 | sed 's/v//'"), topic: versions, emit: versions_blobtoolkit

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1 &&
        (workflow.profile.tokenize(',').intersect(['docker', 'singularity', 'podman', 'apptainer']).size() == 0)) {
        exit 1, "BLOBTOOLKIT_WINDOWSTATS module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    btk pipeline window-stats \\
            --in ${tsv} \\
            $args \\
            --out ${prefix}_window_stats.tsv
    """
}
