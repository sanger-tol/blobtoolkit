process UNCHUNK_BLASTX {
    tag "$meta.id"
    label 'process_single'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the UNCHUNK_BLASTX module. Please use docker or singularity containers."
    }
    container 'genomehubs/blobtoolkit:4.1.2'

    input:
    tuple val(meta), path(raw_proteomes)

    output:
    tuple val(meta), path('*.reference_proteomes.out') , emit: proteomes
    path "versions.yml"                                , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    btk pipeline unchunk-blast \\
        $args \\
        --in ${raw_proteomes} \\
        --out ${prefix}.reference_proteomes.out
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
         blobtoolkit: \$(btk --version | cut -d' ' -f2 | sed 's/v//')
    END_VERSIONS
    """
}
