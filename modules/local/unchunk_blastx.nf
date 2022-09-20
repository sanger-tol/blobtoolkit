process UNCHUNK_BLASTX {
    tag "$meta.id"

    container "genomehubs/blobtoolkit-blobtools:3.3.4"

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
