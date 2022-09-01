//
//
// Run BUSCO for a genome from GOAT and runs diamond_blastp
//
//

nextflow.enable.dsl = 2

include { GOAT_TAXONSEARCH    } from '../../modules/local/goat_taxon_search'
include { BUSCO               } from '../../modules/nf-core/modules/busco/main'
//include { EXTRACT_BUSCO_GENES } from '../../modules/local/extract_busco_genes'
//include { DIAMOND_BLASTP      } from '../../modules/nf-core/modules/diamond/blastp/main'

workflow BUSCO_DIAMOND {
    take:
    //  Tuple [meta, fasta]:
    fasta

    main:

    ch_versions = Channel.empty()

    // this is the string used to name all intermediate and final output files
    name = fasta.map { f -> f.simpleName }
    //meta_id = name.first().toString()

    //
    // Fetch BUSCO lineages for taxon (or taxa)
    //

    GOAT_TAXONSEARCH (
    fasta.map { fa -> [fa[0], params.taxon, params.taxa_file ? file(params.taxa_file) : []] }
    )
    ch_versions = ch_versions.mix(GOAT_TAXONSEARCH.out.versions)

    //
    // Run BUSCO search
    //

    ch_busco_inputs = fasta.join(GOAT_TAXONSEARCH.out.busco_lineages.map { f -> f.readLines().join(',') }) // readLines() transforms all lines to a list

    BUSCO (
    ch_busco_inputs.map { it[0] },
    ch_busco_inputs.map { it[1] },
    params.busco_lineages_path,  // Please pass this option. We don't want to download the lineage data every time.
    [] // No config
    )
    ch_versions = ch_versions.mix(BUSCO.out.versions)

    //
    // Extract BUSCO genes
    //

    //dir = BUSCO.out.busco_dir.map { fileid, path -> path }
    //tables = Channel.fromPath( ["$dir/**/run_archaea_odb10/full_table.tsv", "$dir/**/run_bacteria_odb10/full_table.tsv", "$dir/**/run_eukaryota_odb10/full_table.tsv"] )
    //busco_tables = tables.toList()
    //EXTRACT_BUSCO_GENES (
     //name,
     //busco_tables
    //)
    //ch_versions = ch_versions.mix(EXTRACT_BUSCO_GENES.out.versions)

    //
    // Runs diamond_blastp with the extracted busco genes
    //

    //DIAMOND_BLASTP (
    //[ [ id:name ],  EXTRACT_BUSCO_GENES.out.fasta ],
    //diamonddb,
    //outext,
    //blast_cols
    //)
    //ch_versions = ch_versions.mix(DIAMOND_BLASTP.out.versions)

    emit:

    // diamond_blastp outputs
    //txt      = DIAMOND_BLASTP.out.txt

    // tool versions
    versions = ch_versions
}
