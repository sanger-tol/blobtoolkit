process EXTRACT_NOHIT_FASTA {
    tag "$fasta"

    container "quay.io/biocontainers/seqtk:1.3--h5bf99c6_3"

    input:
    path fasta
    path blastx
    val  evalue

    output:
    path "*.nohit.fasta" , emit: nohit_fasta
    path "versions.yml"  , emit: versions

    script:
    // should check if awk is available within container,
    // if not use separate module (input: blastx, evalue) that uses an image with AWK
    def fname = fasta.simpleName()
    """
    seqtk subseq ${fasta} <(grep '>' ${fasta} | \\
        grep -v -w -f <(awk '{{if($14<{"$evalue"}){{print $1}}}}' ${blastx} | \\
        sort | uniq) \\
        | cut -f1 | sed 's/>//') > "$fname".nohit.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
