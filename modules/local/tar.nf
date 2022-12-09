process TAR {
    tag "$meta.id"

    container "ubuntu:latest"

    input:
    tuple val(meta), path(tbl_a), path(tbl_b), path(tbl_e)
    tuple val(meta), path(seq_a), path(seq_b), path(seq_e)

    output:
    tuple val(meta), path('archaea_odb10/*')   , emit: dir_a
    tuple val(meta), path('bacteria_odb10/*')  , emit: dir_b
    tuple val(meta), path('eukaryota_odb10/*') , emit: dir_e
    path "versions.yml"                        , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def parent_a = "archaea_odb10"
    def parent_b = "bacteria_odb10"
    def parent_e = "eukaryota_odb10"
    """
    mkdir $parent_a
    mkdir $parent_b
    mkdir $parent_e
    gzip --no-name $tbl_a; mv "$tbl_a".gz $parent_a/full_table.tsv.gz
    gzip --no-name $tbl_b; mv "$tbl_b".gz $parent_b/full_table.tsv.gz
    gzip --no-name $tbl_e; mv "$tbl_e".gz $parent_e/full_table.tsv.gz
    tar -zcvf "$seq_a".tar.gz $seq_a; mv "$seq_a".tar.gz $parent_a/busco_sequences.tar.gz
    tar -zcvf "$seq_b".tar.gz $seq_b; mv "$seq_b".tar.gz $parent_b/busco_sequences.tar.gz
    tar -zcvf "$seq_e".tar.gz $seq_e; mv "$seq_e".tar.gz $parent_e/busco_sequences.tar.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tar: \$(tar --version | cut -d' ' -f2 | sed 's/v//')
    END_VERSIONS
    """
}
