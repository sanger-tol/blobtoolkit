process BLOBTOOLKIT_UNCHUNK {
    tag "$meta.id"
    label 'process_single'

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "BLOBTOOLKIT_UNCHUNK module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    container "docker.io/genomehubs/blobtoolkit:4.3.9"

    input:
    tuple val(meta), path(blast_table)

    output:
    tuple val(meta), path("*.out"), emit: blast_out
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${blast_table}"
    """
    btk pipeline unchunk-blast \\
        --in ${blast_table} \\
        --out ${prefix}.out \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: \$(btk --version | cut -d' ' -f2 | sed 's/v//')
    END_VERSIONS
    """
}
