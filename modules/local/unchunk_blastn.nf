process UNCHUNK_BLASTN {
    tag "$meta.id"

    container "genomehubs/blobtoolkit-blobtools:3.3.4"

    input:
    tuple val(meta), path(raw_blastn)

    output:
    tuple val(meta), path('*.blastn') , emit: blastn_out
    path "versions.yml"               , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    btk pipeline unchunk-blast \\
         $args \\
        --in ${raw_blastn} \\
        --out ${prefix}.blastn

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: \$( echo "3.1.0" | sed 's/blobtoolkit //g')
    END_VERSIONS
    """
}
