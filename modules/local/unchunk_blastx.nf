process CHUNK_FASTA_BUSCO {
    tag "$fasta"

    container "genomehubs/blobtoolkit-blobtools"

    input:
    path raw_proteomes
    val max_target_seqs

    output:
    path "reference_proteomes.out" , emit: proteomes
    path "versions.yml"            , emit: versions

    script:
    """
    btk pipeline unchunk-blast \\
        --in ${raw_proteomes} \\
        --count $max_target_seqs \\
        --out reference_proteomes.out 2> unchunk_blastx.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: \$( echo "3.1.0" | sed 's/blobtoolkit //g')
    END_VERSIONS
    """
}
