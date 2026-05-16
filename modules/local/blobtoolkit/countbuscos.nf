process BLOBTOOLKIT_COUNTBUSCOS {
    tag "${meta2.id}"
    label 'process_single'

    container "docker.io/genomehubs/blobtoolkit:4.4.6"

    input:
    tuple val(meta), path(table, stageAs: 'dir??/*')
    tuple val(meta2), path(bed)

    output:
    tuple val(meta2), path("*_buscogenes.tsv"), emit: tsv
    tuple val("${task.process}"), val("blobtoolkit"), eval("btk --version | cut -d' ' -f2 | sed 's/v//'"), topic: versions, emit: versions_blobtoolkit

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1 {
        exit 1, "BLOBTOOLKIT_COUNTBUSCOS module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta2.id}"
    def busco_inputs = (table instanceof List ? table : [table]).collect{ file -> "--in " + file }.join(' ')
    """
    btk pipeline count-busco-genes \\
        $busco_inputs \\
        --mask ${bed} \\
        --out ${prefix}_buscogenes.tsv \\
        ${args}
    """
}
