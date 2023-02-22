//
//
// Generate .yaml config file from meta, then run (1) create, (2) add, and (3) add_summary_to_metadata
//
//

nextflow.enable.dsl = 2

include { GENERATE_CONFIG } from '../../modules/local/generate_config'

workflow BLOBTOOLS {
    take:

    //  Tuple [meta, fasta]:
    fasta

    main:

    ch_versions = Channel.empty()

    //
    // Generate config file
    //
    GENERATE_CONFIG (
    fasta.map { it[0] }
    )
    ch_versions = ch_versions.mix(GENERATE_CONFIG.out.versions)  

    emit: 

    // YAML config file
    config = GENERATE_CONFIG.out.yaml
    
    // tool versions
    versions = ch_versions
}
