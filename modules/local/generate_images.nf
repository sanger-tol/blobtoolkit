process GENERATE_IMAGES {
    tag "$meta.id"
    label 'process_single'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the GENERATE_IMAGES module. Please use docker or singularity containers."
    }
    container 'genomehubs/blobtoolkit:4.0.7'

    input:
    tuple val(meta), path(blobdir)

    output:
    tuple val(meta), path('*.png') , emit: png
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    if( params.use_cov ) {
    """
    blobtools view --view blob --param plotShape=circle --param largeFonts=true --format png --out ${blobdir} "${blobdir}" $args
    blobtools view --view blob --param plotShape=hex --param largeFonts=true --format png --out ${blobdir} "${blobdir}" $args 
    blobtools view --view blob --param plotShape=square --param largeFonts=true --format png --out ${blobdir} "${blobdir}" $args
    blobtools view --view blob --param plotShape=kite --param largeFonts=true --format png --out ${blobdir} "${blobdir}" $args
    blobtools view --view cumulative --param largeFonts=true --format png --out ${blobdir} "${blobdir}" $args
    blobtools view --view snail --param largeFonts=true --format png --out ${blobdir} "${blobdir}" $args
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit-pipeline: \$(blobtoolkit-pipeline --version | cut -d' ' -f2 | sed 's/v//')
    END_VERSIONS
    """
    }
    else {
    """
    blobtools view --view cumulative --param largeFonts=true --format png --out ${blobdir} "${blobdir}" $args
    blobtools view --view snail --param largeFonts=true --format png --out ${blobdir} "${blobdir}" $args
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtoolkit-pipeline: \$(blobtoolkit-pipeline --version | cut -d' ' -f2 | sed 's/v//')
    END_VERSIONS
    """
    }
}
