process BLOBTOOLKIT_CONFIG {
    tag "$meta.id"
    label 'process_single'

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "GENERATE_CONFIG module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    container "docker.io/genomehubs/blobtoolkit:4.3.3"

    input:
    tuple val(meta), val(reads)
    tuple val(meta), val(fasta)

    output:
    tuple val(meta), path("${meta.id}/*.yaml"), emit: yaml
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_reads = reads.collect{"--reads $it"}.join(' ')
    """
    btk pipeline \\
        generate-config \\
        ${prefix} \\
        $args \\
        ${input_reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: \$(btk --version | cut -d' ' -f2 | sed 's/v//')
    END_VERSIONS
    """
}
