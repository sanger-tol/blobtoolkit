process RESTRUCTUREBUSCODIR {
    tag "${meta.id}_${lineage}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ubuntu:20.04'
        : 'docker.io/ubuntu:20.04'}"

    input:
    tuple val(meta), val(lineage), path(batch_summary), path(short_summary_txt), path(short_summary_json), path(full_table), path(missing_busco_list), path(seq_dir)

    output:
    tuple val(meta), path("${lineage}"), emit: clean_busco_dir
    tuple val("${task.process}"), val('restructurebuscodir'), eval('echo 1.0'), topic: versions, emit: versions_restructurebuscodir

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p ${lineage}

    [ -e "${batch_summary}" ]      && ln -s ../${batch_summary}      ${lineage}/short_summary.tsv
    [ -e "${short_summary_txt}" ]  && ln -s ../${short_summary_txt}  ${lineage}/short_summary.txt
    [ -e "${short_summary_json}" ] && ln -s ../${short_summary_json} ${lineage}/short_summary.json
    [ -e "${full_table}" ]         && ln -s ../${full_table}         ${lineage}/
    [ -e "${missing_busco_list}" ] && ln -s ../${missing_busco_list} ${lineage}/

    if [ -e "${seq_dir}/single_copy_busco_sequences" ]
    then
        tar czf ${lineage}/single_copy_busco_sequences.tar.gz -C ${seq_dir} single_copy_busco_sequences
        tar czf ${lineage}/multi_copy_busco_sequences.tar.gz  -C ${seq_dir} multi_copy_busco_sequences
        tar czf ${lineage}/fragmented_busco_sequences.tar.gz  -C ${seq_dir} fragmented_busco_sequences
    else
        # Busco was run in --tar mode, the sequences are already compressed
        ln -s ../${seq_dir}/single_copy_busco_sequences.tar.gz ${lineage}/
        ln -s ../${seq_dir}/multi_copy_busco_sequences.tar.gz ${lineage}/
        ln -s ../${seq_dir}/fragmented_busco_sequences.tar.gz ${lineage}/
    fi
    """

    stub:
    """
    mkdir -p ${lineage}
    touch ${lineage}/short_summary.tsv
    touch ${lineage}/short_summary.txt
    touch ${lineage}/short_summary.json
    touch ${lineage}/full_table.tsv
    touch ${lineage}/missing_busco_list.tsv
    tar -czf ${lineage}/single_copy_busco_sequences.tar.gz -T /dev/null
    tar -czf ${lineage}/multi_copy_busco_sequences.tar.gz  -T /dev/null
    tar -czf ${lineage}/fragmented_busco_sequences.tar.gz  -T /dev/null
    """
}
