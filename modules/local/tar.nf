process TAR {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "conda-forge::coreutils=9.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    tuple val(meta), path(tbl_a), path(tbl_b), path(tbl_e)
    tuple val(meta), path(seq_a), path(seq_b), path(seq_e)

    output:
    tuple val(meta), path('archaea_odb10'), path('bacteria_odb10'), path('eukaryota_odb10'), emit: dir_abe
    path "versions.yml"                                                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def parent_a = "archaea_odb10"
    def parent_b = "bacteria_odb10"
    def parent_e = "eukaryota_odb10"
    """
    mkdir $parent_a
    mkdir $parent_b
    mkdir $parent_e
    gzip --no-name --force $tbl_a; mv "$tbl_a".gz $parent_a/full_table.tsv.gz
    gzip --no-name --force $tbl_b; mv "$tbl_b".gz $parent_b/full_table.tsv.gz
    gzip --no-name --force $tbl_e; mv "$tbl_e".gz $parent_e/full_table.tsv.gz
    tar -zcvf "$seq_a".tar.gz $seq_a; mv "$seq_a".tar.gz $parent_a/busco_sequences.tar.gz
    tar -zcvf "$seq_b".tar.gz $seq_b; mv "$seq_b".tar.gz $parent_b/busco_sequences.tar.gz
    tar -zcvf "$seq_e".tar.gz $seq_e; mv "$seq_e".tar.gz $parent_e/busco_sequences.tar.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tar: \$(echo \$(tar --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
