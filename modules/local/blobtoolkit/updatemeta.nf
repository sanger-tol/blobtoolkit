process BLOBTOOLKIT_UPDATEMETA {
    tag "$meta.id"
    label 'process_single'

    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "BLOBTOOLKIT_UPDATEMETA module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    container "docker.io/genomehubs/blobtoolkit:4.3.9"

    input:
    tuple val(meta), path(input)
    val reads
    path versions
    // The following are passed as "val" because we just want to know the full paths. No staging necessary
    val blastp
    val blastx
    val blastn
    val taxdump

    output:
    tuple val(meta), path("*.json"), emit: json
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def input_reads = reads.collect{"--read_id ${it[0].id} --read_type ${it[0].datatype} --read_path ${it[1]}"}.join(' ')
    """
    update_versions.py \\
        ${args} \\
        --meta_in ${input}/meta.json \\
        --software ${versions} \\
        --blastp ${blastp} \\
        --blastx ${blastx} \\
        --blastn ${blastn} \\
        --taxdump ${taxdump} \\
        $input_reads \\
        --meta_out ${prefix}.meta.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        update_versions.py: \$(update_versions.py --version | cut -d' ' -f2)
    END_VERSIONS
    """

}
