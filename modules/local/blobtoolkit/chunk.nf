process BLOBTOOLKIT_CHUNK {
    tag "$meta.id"
    label 'process_single'

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "BLOBTOOLKIT_CHUNK module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    container "docker.io/genomehubs/blobtoolkit:4.3.13"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(busco_table)

    output:
    tuple val(meta), path("*.chunks.fasta"), emit: chunks
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def busco = busco_table ? "--busco ${busco_table}" : "--busco None"
    """
    btk pipeline chunk-fasta \\
        --in ${fasta} \\
        ${busco} \\
        --out ${prefix}.chunks.fasta \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: \$(btk --version | cut -d' ' -f2 | sed 's/v//')
    END_VERSIONS
    """
}
