process BLOBTOOLKIT_COUNTBUSCOS {
    tag "$meta.id"
    label 'process_single'

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "BLOBTOOLKIT_COUNTBUSCOS module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    container "genomehubs/blobtoolkit:4.1.2"

    input:
    tuple val(meta), path(table, stageAs: 'dir??/*')
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*_buscogenes.tsv"), emit: tsv
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def busco_inputs = table.collect{"--in $it"}.join(' ')
    """
    btk pipeline count-busco-genes \\
        $busco_inputs \\
        --mask ${bed} \\
        --out ${prefix}_buscogenes.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: \$(btk --version | cut -d' ' -f2 | sed 's/v//')
    END_VERSIONS
    """
}
