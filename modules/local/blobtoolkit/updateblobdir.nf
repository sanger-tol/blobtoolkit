process BLOBTOOLKIT_UPDATEBLOBDIR {
    tag "$meta.id"
    label 'process_medium'

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "BLOBTOOLKIT_BLOBDIR module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    container "genomehubs/blobtoolkit:4.1.5"

    input:
    tuple val(meta), path(input)
    tuple val(meta1), path(blastx)
    path(taxdump)

    output:
    tuple val(meta), path(prefix), emit: blobdir
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def hits_blastx = blastx ? "--hits ${blastx}" : ""
    """
    blobtools replace \\
        --taxdump ${taxdump} \\
        --taxrule bestdistorder=buscoregions \\
        ${hits_blastx} \\
        --bedtsvdir windowstats \\
        --meta ${yaml} \\
        --threads ${task.cpus} \\
        $args \\
        ${input}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: \$(btk --version | cut -d' ' -f2 | sed 's/v//')
    END_VERSIONS
    """
}
