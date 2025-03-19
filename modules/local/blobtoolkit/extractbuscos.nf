process BLOBTOOLKIT_EXTRACTBUSCOS {
    tag "$meta.id"
    label 'process_single'

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "BLOBTOOLKIT_EXTRACTBUSCOS module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    container "docker.io/genomehubs/blobtoolkit:4.4.4"

    input:
    tuple val(meta), path(fasta)
    path seq, stageAs: "lineage??/*"

    output:
    tuple val(meta), path("*_buscogenes.fasta"), emit: genes
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def seq_args = (seq instanceof List ? seq : [seq]).collect { "--busco " + it } .join(' ')
    """
    btk pipeline extract-busco-genes \\
        $seq_args \\
        --out ${prefix}_buscogenes.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: \$(btk --version | cut -d' ' -f2 | sed 's/v//')
    END_VERSIONS
    """
}
