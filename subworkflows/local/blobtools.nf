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
    // Generate or read config file: a YAML file or a GCA accesion should be provided
    //
    if ( params.accesion && !params.yaml){
      GENERATE_CONFIG (
      fasta.map { fa -> [fa[0], "${params.accesion}"] }
      )
      ch_versions = ch_versions.mix(GENERATE_CONFIG.out.versions)
      config_file = GENERATE_CONFIG.out.yaml
    }
    if ( params.yaml && !params.accesion){
      config_file = Channel.fromPath(params.yaml)
    }
    if ( (!params.accesion && !params.yaml) || (params.accesion && params.yaml) ){
      exit 1, 'Input not specified. Please include either a YAML file for draft genome or GCA accesion for published genome'
    }
   
    emit:

    // YAML config file
    config_file
    
    // tool versions
    versions = ch_versions
}
