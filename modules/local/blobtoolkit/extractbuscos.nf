process BLOBTOOLKIT_EXTRACTBUSCOS {
    tag "$meta.id"
    label 'process_single'

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "BLOBTOOLKIT_EXTRACTBUSCOS module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    container "genomehubs/blobtoolkit:4.1.2"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta1), path(table1, stageAs: "lineage1/*"), path(seq1, stageAs: "lineage1/*")
    tuple val(meta2), path(table2, stageAs: "lineage2/*"), path(seq2, stageAs: "lineage2/*")
    tuple val(meta3), path(table3, stageAs: "lineage3/*"), path(seq3, stageAs: "lineage3/*")

    output:
    tuple val(meta), path("*_buscogenes.fasta"), emit: genes
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    btk pipeline extract-busco-genes \\
        --busco $table1 \\
        --busco $table2 \\
        --busco $table3 \\
        --out ${prefix}_buscogenes.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: \$(btk --version | cut -d' ' -f2 | sed 's/v//')
    END_VERSIONS
    """
}
