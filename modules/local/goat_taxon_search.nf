process GOAT_TAXONSEARCH {
    tag "$busco_lineages"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::goat=0.2.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/goat:0.2.0--h92d785c_0':
        'quay.io/biocontainers/goat:0.2.0--h92d785c_0' }"

    input:
    val taxon
    path taxa_file

    output:
    path "*.tsv" 	    , emit: taxonsearch
    file "*.txt" 	    , emit: busco_lineages
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ''
    input = taxa_file ? "-f ${taxa_file}" : "-t ${taxon}"
    if (!taxon && !taxa_file) error "No input. Valid input: single taxon identifier or a .txt file with identifiers"
    if (taxon && taxa_file ) error "Only one input is required: a single taxon identifier or a .txt file with identifiers"
    // ${prefix}.txt contains the list of BUSCO (odb10) lineages, one lineage per line without empty lines
    """
    goat-cli taxon search \\
        $args \\
        $input > ${prefix}.tsv

    cat ${prefix}.tsv | cut -f5 | sed '1d' | grep . > ${prefix}.txt
    echo "bacteria_odb10" >> ${prefix}.txt
    echo "archaea_odb10" >> ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        goat: \$(goat-cli --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
