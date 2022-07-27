``//
// Run BUSCO for a genome from GOAT and runs diamond_blastp
//

include { GOAT_TAXONSEARCH    } from '../../modules/local/goat_taxon_search'
include { EXTRACT_BUSCO_GENES } from '../../modules/local/extract_busco_genes'
include { BUSCO               } from '../../modules/nf-core/modules/busco/main'
include { DIAMOND_BLASTP      } from '../../modules/nf-core/modules/diamond/blastp/main'

workflow BUSCO_DIAMOND {
    take:

    // GOAT_TAXONSEARCH input
    taxon // value: single binomial name or NCBI taxonomy ID or '' if taxa_file is provided
    taxa_file // file containing a taxon ID per line or empty list if taxon is provided

    // BUSCO input
    genome_fasta //  path to genome fasta file (s)

    // diamond_blastp input
    diamond_db // path to diamond database
    out_ext // module parameter
    blast_columns // module parameter

    main:

    ch_versions = Channel.empty()

    //
    // Fetch BUSCO lineages for taxon (or taxa)
    //
    GOAT_TAXONSEARCH (
        [ taxon, taxa_file ] // input: taxon, no taxa_file or no taxon, taxa_file
    )
    ch_versions = ch_versions.mix(GOAT_TAXONSEARCH.out.versions.first())

    //
    // Run BUSCO search
    //
    BUSCO (
    [ [ id:'blobtoolkit' ],  file(genome_fasta) ],
    GOAT_TAXONSEARCH.out.busco_lineages.readLines(), // readLines() transforms all lines to a list
    [], // Download busco lineage
    []  // No config file, in this case it might be required
    )
    ch_versions = ch_versions.mix(BUSCO.out.versions)

    //
    // Extract BUSCO genes
    //
    EXTRACT_BUSCO_GENES (
     busco_table // should be copied from busco_dir output from BUSCO and transformed.
    )
    ch_versions = ch_versions.mix(EXTRACT_BUSCO_GENES.out.versions)

    //
    // Runs diamond_blastp with the extracted busco genes
    //
    DIAMOND_BLASTP (
    [ [ id:'blobtoolkit' ],  EXTRACT_BUSCO_GENES.out.fasta ],
    diamond_db,
    out_ext,
    blast_columns
    )
    ch_versions = ch_versions.mix(DIAMOND_BLASTP.out.versions)

    emit:

    // check if there are no missing outputs
    // or outputs that will be not used later

    // diamond_blastp outputs
    blast    = DIAMOND_BLASTP.out.blast
    xml      = DIAMOND_BLASTP.out.xml
    txt      = DIAMOND_BLASTP.out.txt
    daa      = DIAMOND_BLASTP.out.daa
    sam      = DIAMOND_BLASTP.out.sam
    tsv      = DIAMOND_BLASTP.out.tsv
    paf      = DIAMOND_BLASTP.out.paf

    versions = ch_versions
}
