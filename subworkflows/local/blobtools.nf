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

    fasta               //  Tuple [meta, path/to/fasta]
    windowstats_tsv     //  Tuple [meta, tuple[ path/to/window_stats.0.01, path/to/window_stats.01, path/to/window_stats.100000, path/to/window_stats.1000000 ]]
    busco_table         //  Tuple [meta, path/to/busco_first_table]
    blastp              //  Tuple [meta, path/to/blastp]
    blobdir_name        //  Val(genome_meta)
    ncbi_taxdump        //  Path(ncbi_taxdump)

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
   ch_versions = ch_versions.mix(ADD_SUMMARY_TO_METADATA.out.versions)

    //  
    // Create Blobdir data structure
    //
    CREATE_BLOBDIR (windowstats_tsv, busco_table, blastp, ADD_SUMMARY_TO_METADATA.out.yaml, ncbi_taxdump, blobdir_name)
    ch_versions = ch_versions.mix(CREATE_BLOBDIR.out.versions)

    emit:

    // YAML config file
    config_yaml = ADD_SUMMARY_TO_METADATA.out.yaml

    // Blobdir directory
    blobdir = CREATE_BLOBDIR.out.blobdir

    // tool versions
    versions = ch_versions
}
