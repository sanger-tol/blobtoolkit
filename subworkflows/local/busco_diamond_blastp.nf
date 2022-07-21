//
// Run BUSCO for a genome from GOAT and runs diamond_blastp
//

include { GOAT_TAXONSEARCH    } from '../../modules/local/goat_taxon_search'
include { EXTRACT_BUSCO_GENES } from '../../modules/local/extract_busco_genes'
include { BUSCO               } from '../../modules/nf-core/modules/busco/main'
include { DIAMOND_BLASTP      } from '../../modules/nf-core/modules/diamond/blastp/main'

workflow BUSCO_DIAMOND {
    take:
    taxon // value: binomial name or NCBI taxonomy ID
    taxa_file // file: /path/to/*taxa_file.txt
    busco_lineages // file: /path/to/busco_lineages.txt
    fasta // file: /path/to/output_busco_genes.fasta

    main:
    // main workflow

    emit:
    // write outputs

}
