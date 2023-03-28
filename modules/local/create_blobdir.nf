process CREATE_BLOBDIR {
    tag "$meta.id"
    label 'process_single'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the CREATE_BLOBDIR module. Please use docker or singularity containers."
    }
    container "genomehubs/blobtoolkit-blobtools:3.3.4"

    input:
    tuple val(meta), path(window, stageAs: 'windowstats/*')
    tuple val(meta), path(busco)
    tuple val(meta), path(blastp)
    tuple val(meta), path(yaml)
    path(taxdump)
    val(genome_accession)

    output:
    tuple val(meta), path('**/*meta.json')              , emit: json
    tuple val(meta), path('**/*identifiers.json')       , emit: identifiers
    tuple val(meta), path('**/*buscogenes_phylum.json') , emit: buscogenes_phylum
    tuple val(meta), path('**/*cov.json')               , emit: coverage
    tuple val(meta), path('**/*busco.json')             , emit: busco
    tuple val(meta), path('**/*.json')                  , emit: blobdir
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    blobtools replace \\
        --bedtsvdir windowstats \\
        --meta ${yaml} \\
        --taxdump ${taxdump} \\
        --taxrule buscogenes \\
        --busco ${busco} \\
        --hits ${blastp} \\
        --threads ${task.cpus} \\
        $args \\
        ${genome_accession}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit: \$(btk --version | cut -d' ' -f2 | sed 's/v//')
    END_VERSIONS
    """
}
