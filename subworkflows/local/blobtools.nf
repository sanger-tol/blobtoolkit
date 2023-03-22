//
//
// Generate .yaml config file from meta, then run (1) create, (2) add, and (3) add_summary_to_metadata
//
//

nextflow.enable.dsl = 2

include { GENERATE_CONFIG         } from '../../modules/local/generate_config'
include { ADD_SUMMARY_TO_METADATA } from '../../modules/local/add_summary_to_metadata'
include { CREATE_BLOBDIR          } from '../../modules/local/create_blobdir'

workflow BLOBTOOLS {
    take:

    //  Tuple [meta, fasta]:
    fasta
    windowstats_tsv
    busco_table
    blastp
    ncbi_taxdump
    blobdir_name

    main:

    ch_versions = Channel.empty()

    //
    // Generate or read config file: a YAML file or a GCA accession should be provided
    //
    if ( params.accession && !params.yaml){
      GENERATE_CONFIG (
      fasta.map { fa -> [fa[0], "${params.accession}"] }
      )
      ch_versions = ch_versions.mix(GENERATE_CONFIG.out.versions)
      config_file = GENERATE_CONFIG.out.yaml
    }
    if ( params.yaml && !params.accession){
      config_file = fasta.map { fa -> [fa[0], "${params.yaml}"] }
    }
    if ( (!params.accession && !params.yaml) || (params.accession && params.yaml) ){
      exit 1, 'Input not specified. Please include either a YAML file for draft genome or GCA accession for published genome'
    }

    //
    // Add summary to metadata
    //
    ADD_SUMMARY_TO_METADATA (
      config_file
    )

    //windowstats_tsv.view()
    //window_dir = GetDirectory(windowstats_tsv)
    //window_dir.view()

    //tsv = TextToTsv.(blastp)

    //  
    // Create Blobdir data structure
    //
    CREATE_BLOBDIR (windowstats_tsv, busco_table, blastp, ncbi_taxdump, ADD_SUMMARY_TO_METADATA.out.yaml, blobdir_name)
   
    emit:

    // YAML config file
    config_yaml = ADD_SUMMARY_TO_METADATA.out.yaml
    json = CREATE_BLOBDIR.out.json
    
    // tool versions
    versions = ch_versions
}
