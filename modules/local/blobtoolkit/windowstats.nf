process BLOBTOOLKIT_WINDOWSTATS {
    tag "$meta.id"

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "BLOBTOOLKIT_WINDOWSTATS module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    container "docker.io/genomehubs/blobtoolkit:develop"

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path('*_window_stats*.tsv') , emit: tsv
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    btk pipeline window-stats \\
            --in ${tsv} \\
            $args \\
            --out ${prefix}_window_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: \$(btk --version | cut -d' ' -f2 | sed 's/v//')
    END_VERSIONS
    """
}
