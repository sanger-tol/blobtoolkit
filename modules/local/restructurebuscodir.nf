process RESTRUCTUREBUSCODIR {
    tag "${meta.id}_${lineage}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), val(lineage), path(batch_summary), path(short_summary_txt), path(short_summary_json), path(full_table), path(missing_busco_list), path(single_copy_busco_sequences), path(multi_copy_busco_sequences), path(fragmented_busco_sequences), path(hmmer_output)

    output:
    tuple val(meta), path("${lineage}"), emit: clean_busco_dir
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${lineage}

    [ -n "${batch_summary}" ] && cp --dereference ${batch_summary} ${lineage}/short_summary.tsv
    [ -n "${short_summary_txt}" ] && cp --dereference ${short_summary_txt} ${lineage}/short_summary.txt
    [ -n "${short_summary_json}" ] && cp --dereference ${short_summary_json} ${lineage}/short_summary.json

    [ -e ${full_table} ] && cp ${full_table} ${lineage}/
    [ -e ${missing_busco_list} ] && cp ${missing_busco_list} ${lineage}/

    tar czf ${lineage}/single_copy_busco_sequences.tar.gz -C \$(dirname ${single_copy_busco_sequences}) \$(basename ${single_copy_busco_sequences})
    tar czf ${lineage}/multi_copy_busco_sequences.tar.gz -C \$(dirname ${multi_copy_busco_sequences}) \$(basename ${multi_copy_busco_sequences})
    tar czf ${lineage}/fragmented_busco_sequences.tar.gz -C \$(dirname ${fragmented_busco_sequences}) \$(basename ${fragmented_busco_sequences})
    tar czf ${lineage}/hmmer_output.tar.gz --exclude=.checkpoint -C \$(dirname ${hmmer_output}) \$(basename ${hmmer_output})

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tar: \$(tar --version | awk 'NR==1 {print \$3}')
    END_VERSIONS
    """
}
