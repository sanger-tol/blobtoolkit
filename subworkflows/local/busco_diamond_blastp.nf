//
//
// Run BUSCO for a genome from GOAT and runs diamond_blastp
//
//

nextflow.enable.dsl = 2

include { GOAT_TAXONSEARCH    } from '../../modules/local/goat_taxon_search'
include { BUSCO               } from '../../modules/nf-core/busco/main'
include { EXTRACT_BUSCO_GENES } from '../../modules/local/extract_busco_genes'
include { DIAMOND_BLASTP      } from '../../modules/nf-core/diamond/blastp/main'

workflow BUSCO_DIAMOND {
    take:

    //  Tuple [meta, fasta]:
    fasta

    main:

    ch_versions = Channel.empty()


    //
    // Fetch BUSCO lineages for taxon (or taxa)
    //

    GOAT_TAXONSEARCH (
    fasta.map { fa -> [fa[0], "${params.taxon}", "${params.taxa_file}" ? file("${params.taxa_file}") : []] }
    )
    ch_versions = ch_versions.mix(GOAT_TAXONSEARCH.out.versions)

    //
    // Run BUSCO search
    //

    // Make a new channel in which each output file is replaced by its content, with meta propagated
    ch_lineages = GOAT_TAXONSEARCH.out.busco_lineages.flatMap { it[1].readLines().collect { line -> [it[0], line] } }
    // Cross-product of both channels, using meta as the key
    ch_busco_inputs = fasta.combine(ch_lineages, by: 0)

    BUSCO (
    ch_busco_inputs.map { [it[0], it[1]] },
    ch_busco_inputs.map { it[2] },
    "${params.busco_lineages_path}",  // Please pass this option. We don't want to download the lineage data every time.
    [] // No config
    )
    ch_versions = ch_versions.mix(BUSCO.out.versions)

    //
    // Extract BUSCO genes
    //

    // gets file name which is a folder name in paths to busco full tables
    fasta_filename = fasta.map { meta,fa -> [meta, fa.name] }

    // channel: emits paths to busco results for each lineage
    dir = BUSCO.out.busco_dir.combine(fasta_filename, by:0)

    // filter paths to busco full tables for archaea, bacteria and eukaryota
    dir_a = dir.filter { "$it" =~ /archaea_odb10/ }.map { meta,a,fname -> [meta, "$a/$fname/run_archaea_odb10/full_table.tsv"] }.collect()
    dir_b = dir.filter { "$it" =~ /bacteria_odb10/ }.map { meta,b,fname -> [meta, "$b/$fname/run_bacteria_odb10/full_table.tsv"] }.collect()
    dir_e = dir.filter { "$it" =~ /eukaryota_odb10/ }.map { meta,e,fname -> [meta, "$e/$fname/run_eukaryota_odb10/full_table.tsv"] }.collect()

    // create copies of full tables with a lineage identifier, avoids file name collision
    dir_a = dir_a.map { meta,table -> [meta, table.copyTo("$table.parent/archaea_odb10_full_table.tsv")] }.collect()
    dir_b = dir_b.map { meta,table -> [meta, table.copyTo("$table.parent/bacteria_odb10_full_table.tsv")] }.collect()
    dir_e = dir_e.map { meta,table -> [meta, table.copyTo("$table.parent/eukaryota_odb10_full_table.tsv")] }.collect()

    // combine all three channels into a single channel: tuple( meta, a, b, e )
    dir_ab = dir_a.combine(dir_b, by:0)
    dir_abe = dir_ab.combine(dir_e, by:0)

    EXTRACT_BUSCO_GENES (
    dir_abe
    )
    ch_versions = ch_versions.mix(EXTRACT_BUSCO_GENES.out.versions)

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
