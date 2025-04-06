process BLOBTOOLKIT_CREATEBLOBDIR {
    tag "$meta.id"
    label 'process_medium'

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "BLOBTOOLKIT_BLOBDIR module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    container "docker.io/genomehubs/blobtoolkit:4.4.5"

    input:
    tuple val(meta), path(window, stageAs: 'windowstats/*')
    tuple val(meta1), path(busco, stageAs: 'lineage??/*')
    tuple val(meta2), path(blastp)
    tuple val(meta3), path(yaml)
    path(taxdump, stageAs: 'taxdump/taxdump.json')

    output:
    tuple val(meta), path(prefix), emit: blobdir
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def busco_args = (busco instanceof List ? busco : [busco]).collect { "--busco " + it } .join(' ')
    def hits_blastp = blastp ? "--hits ${blastp}" : ""
    """
    blobtools replace \\
        --bedtsvdir windowstats \\
        --meta ${yaml} \\
        --taxdump \$(dirname ${taxdump}) \\
        --taxrule buscogenes \\
        ${busco_args} \\
        ${hits_blastp} \\
        --threads ${task.cpus} \\
        $args \\
        ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: \$(btk --version | cut -d' ' -f2 | sed 's/v//')
    END_VERSIONS
    """
}
