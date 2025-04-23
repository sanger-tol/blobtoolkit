process BLAST_BLASTN {
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.15.0--pl5321h6f7f691_1':
        'biocontainers/blast:2.15.0--pl5321h6f7f691_1' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(db)
    val taxid

    output:
    tuple val(meta), path('*.txt'), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getExtension() == "gz" ? true : false
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta
    def exclude_taxon = taxid ? "-negative_taxids ${taxid}" : ''
    def command_epilog = taxid ? "|| true" : ''

    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    DB=`find -L ./ -name "*.nal" | sed 's/\\.nal\$//'`
    if [ -z "\$DB" ]; then
        DB=`find -L ./ -name "*.nin" | sed 's/\\.nin\$//'`
    fi
    echo Using \$DB

    if [ -n "${taxid}" ]; then
        # Symlink the tax* files (needed for -taxid options to work)
        for file in taxdb.btd taxdb.bti taxonomy4blast.sqlite3; do
            if [ ! -f ${db}/\$file ]; then
                echo "Error: \$file not found in ${db}"
                exit 1
            fi
            ln -s ${db}/\$file .
        done
    fi

    timeout 11.9h blastn \\
      -num_threads ${task.cpus} \\
      -db \$DB \\
      -query ${fasta_name} \\
      ${exclude_taxon} \\
      ${args} \\
      -out ${prefix}.txt \\
      2> >( tee "${prefix}.error.log" >&2 ) $command_epilog

    # Fallback if blastn fails or times out — make sure output exists
    if [ ! -s "${prefix}.txt" ]; then
    echo "blastn failed or timed out — creating empty output"
    touch "${prefix}.txt"
    fi

    if [[ -s "${prefix}.error.log" ]]
    then
        grep -qF 'BLAST Database error: Taxonomy ID(s) not found.Taxonomy ID(s) not found' "${prefix}.error.log"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
