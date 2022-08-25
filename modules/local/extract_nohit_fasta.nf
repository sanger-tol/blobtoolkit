process EXTRACT_NOHIT_FASTA {
    tag "$fasta"

    conda (params.enable_conda ? "bioconda::seqtk=1.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.3--h5bf99c6_3' :
        'quay.io/biocontainers/seqtk:1.3--h5bf99c6_3' }"

    input:
    path fasta
    path nohit_list

    output:
    path "*.nohit.fasta" , emit: nohit_fasta
    path "versions.yml"  , emit: versions

    script:
    def fname = fasta.simpleName()
    """
    seqtk subseq \\
        ${fasta} \\
        ${nohit_list} > "$fname".nohit.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
