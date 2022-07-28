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
    // Value: single binomial name or NCBI taxonomy ID or '' if taxa_file is provided
    taxon = ""
    // File containing a taxon ID per line or empty list if taxon is provided
    taxa_file = []

    // BUSCO input
    //  Path to genome fasta file:
    genome_fasta = []
    // Path to busco lineages - downloads if not set
    lineages_path = []
    // BUSCO configuration file
    busco_config = []

    // diamond_blastp input
    // Directory containing the protein blast database:
    diamonddb = []
    // Specify the type of output file to be generated, `txt` corresponds to to BLAST tabular format:
    outext = "txt"
    // Space separated list of DIAMOND tabular BLAST output keywords:
    blast_cols = "qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

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
    [ [ id:'blobtoolkit' ], genome_fasta ],
    GOAT_TAXONSEARCH.out.busco_lineages.readLines(), // readLines() transforms all lines to a list
    lineages_path,
    busco_config
    )
    ch_versions = ch_versions.mix(BUSCO.out.versions)

    //
    // Extract BUSCO genes
    //
    EXTRACT_BUSCO_GENES (
     BUSCO.out.busco_dir
    )
    ch_versions = ch_versions.mix(EXTRACT_BUSCO_GENES.out.versions)

    //
    // Runs diamond_blastp with the extracted busco genes
    //
    DIAMOND_BLASTP (
    [ [ id:'blobtoolkit' ],  EXTRACT_BUSCO_GENES.out.fasta ],
    diamonddb,
    outext,
    blast_cols
    )
    ch_versions = ch_versions.mix(DIAMOND_BLASTP.out.versions)

    emit:

    // diamond_blastp outputs
    txt      = DIAMOND_BLASTP.out.txt
    // tool versions
    versions = ch_versions
}
