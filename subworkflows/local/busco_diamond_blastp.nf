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
    //  Path to genome fasta file:
    fasta

    // GOAT_TAXONSEARCH input
    /// Value: single binomial name or NCBI taxonomy ID or '' if taxa_file is provided
    taxon
    /// File containing a taxon ID per line or empty list if taxon is provided
    taxa_file

    // DIAMOND_BLASTP input
    /// Directory containing the protein blast database:
    //diamonddb
    /// Specify the type of output file to be generated, `txt` corresponds to to BLAST tabular format:
    //outext
    /// Space separated list of DIAMOND tabular BLAST output keywords:
    /// "qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    //blast_cols

    main:

    ch_versions = Channel.empty()

    // this is the string used to name all intermediate and final output files
    name = fasta.map { f -> f.simpleName }


    //
    // Fetch BUSCO lineages for taxon (or taxa)
    //

    GOAT_TAXONSEARCH (
        name,
        taxon,
        taxa_file // input: taxon, no taxa_file or no taxon, taxa_file
    )
    ch_versions = ch_versions.mix(GOAT_TAXONSEARCH.out.versions)

    //
    // Run BUSCO search
    //

    lineages_list = GOAT_TAXONSEARCH.out.busco_lineages.map { f -> f.readLines() } // readLines() transforms all lines to a list
    lineages = Channel.of( lineages_list )

    BUSCO (
    fasta.map { fa -> [ [id: fa.baseName ], fa ] }, // Add meta data using the file's basename as id
    "${lineages.join()}",
    [], // Download busco lineage
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
