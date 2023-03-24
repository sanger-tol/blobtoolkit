process GET_WINDOW_STATS {
    tag "$meta.id"
    label 'process_single'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the GET_WINDOW_STATS module. Please use docker or singularity containers."
    }
    container "genomehubs/blobtoolkit-blobtools:3.5.4"

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
