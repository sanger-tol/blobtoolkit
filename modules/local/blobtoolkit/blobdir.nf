process BLOBTOOLKIT_BLOBDIR {
    tag "$meta.id"
    label 'process_medium'

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "BLOBTOOLKIT_BLOBDIR module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    container "genomehubs/blobtoolkit:4.1.5"

    input:
    tuple val(meta), path(window, stageAs: 'windowstats/*')
    tuple val(meta1), path(busco)
    tuple val(meta2), path(blastp)
    tuple val(meta3), path(blastx)
    tuple val(meta4), path(blastn)
    tuple val(meta5), path(yaml)
    path(taxdump)

    output:
    tuple val(meta), path(prefix), emit: blobdir
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def hits_blastp = blastp ? "--hits ${blastp}" : ""
    def hits_blastx = blastx ? "--hits ${blastx}" : ""
    def hits_blastn = blastn ? "--hits ${blastn}" : ""
    """
    blobtools replace \\
        --bedtsvdir windowstats \\
        --meta ${yaml} \\
        --taxdump ${taxdump} \\
        --taxrule bestdistorder=buscoregions \\
        --busco ${busco} \\
        ${hits_blastp} \\
        ${hits_blastx} \\
        ${hits_blastn} \\
        --threads ${task.cpus} \\
        $args \\
        ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: \$(btk --version | cut -d' ' -f2 | sed 's/v//')
    END_VERSIONS
    """
}
