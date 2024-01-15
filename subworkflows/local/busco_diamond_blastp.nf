//
// Run BUSCO for a genome from GOAT and runs diamond_blastp
//

include { GOAT_TAXONSEARCH          } from '../../modules/nf-core/goat/taxonsearch/main'
include { BUSCO                     } from '../../modules/nf-core/busco/main'
include { BLOBTOOLKIT_EXTRACTBUSCOS } from '../../modules/local/blobtoolkit/extractbuscos'
include { DIAMOND_BLASTP            } from '../../modules/nf-core/diamond/blastp/main'


workflow BUSCO_DIAMOND {
    take:
    fasta        // channel: [ val(meta), path(fasta) ]
    taxon_taxa   // channel: [ val(meta, val(taxon), path(taxa) ]
    busco_db     // channel: path(busco_db)
    blastp       // channel: path(blastp_db)
    outext       // channel: val(out_format)
    cols         // channel: val(column_names)


    main:
    ch_versions = Channel.empty()


    //
    // Fetch BUSCO lineages for taxon (or taxa)
    //
    GOAT_TAXONSEARCH ( taxon_taxa )
    ch_versions = ch_versions.mix ( GOAT_TAXONSEARCH.out.versions.first() )
    

    //
    // Get NCBI species ID
    //
    GOAT_TAXONSEARCH.out.taxonsearch
    | map { meta, csv -> csv.splitCsv(header:true, sep:'\t', strip:true) }
    | map { row -> [ row.taxon_rank, row.taxon_id ] }
    | transpose()
    | filter { rank,id -> rank =~ /species/ }
    | map { rank, id -> id}
    | set { ch_taxid }


    //
    // Run BUSCO search
    //
    GOAT_TAXONSEARCH.out.taxonsearch
    | map { meta, csv -> csv.splitCsv(header:true, sep:'\t', strip:true) }
    | map { row -> row.odb10_lineage.findAll { it != "" } }
    | set { ch_ancestral_lineages }


    // Add the basal lineages to the list (excluding duplicates)
    basal_lineages = [ "eukaryota_odb10", "bacteria_odb10", "archaea_odb10" ]
    ch_ancestral_lineages
    | map { lineages -> (lineages + basal_lineages).unique() }
    | flatten ()
    | set { ch_lineages }

    BUSCO ( fasta, "genome", ch_lineages, busco_db.collect().ifEmpty([]), [] )
    ch_versions = ch_versions.mix ( BUSCO.out.versions.first() )


    //
    // Select input for BLOBTOOLKIT_EXTRACTBUSCOS
    //
    BUSCO.out.seq_dir
    | filter { meta, seq -> basal_lineages.contains(seq.parent.baseName.minus("run_")) }
    | groupTuple()
    | set { ch_basal_buscos }


    // Extract BUSCO genes from the basal lineages
    BLOBTOOLKIT_EXTRACTBUSCOS ( fasta, ch_basal_buscos )
    ch_versions = ch_versions.mix ( BLOBTOOLKIT_EXTRACTBUSCOS.out.versions.first() )


    //
    // Align BUSCO genes against the BLASTp database
    //    
    BLOBTOOLKIT_EXTRACTBUSCOS.out.genes
    | filter { it[1].size() > 140 }
    | set { ch_busco_genes }

    DIAMOND_BLASTP ( ch_busco_genes, blastp, outext, cols )
    ch_versions = ch_versions.mix ( DIAMOND_BLASTP.out.versions.first() )


    // Index the lineages in the taxonomic order
    def lineage_index = 0
    ch_lineages
    | map { lineage -> [lineage, lineage_index++] }
    | set { ch_ordered_lineages }


    // Order BUSCO results accoring to ch_lineages
    BUSCO.out.full_table
    | map { meta, table -> [table.parent.baseName.minus("run_"), meta, table] }
    | join ( ch_ordered_lineages )
    | map { lineage, meta, table, index -> [meta, table, index] }
    | groupTuple()
    | map { meta, tables, indexes -> [ meta, tables.withIndex().sort { a, b -> indexes[a[1]] <=> indexes[b[1]] } . collect { table, i -> table } ] }
    | set { ch_indexed_buscos }


    // Select BUSCO results for taxonomically closest database
    ch_indexed_buscos
    | map { meta, tables -> [meta, tables[0]] }
    | set { ch_first_table }


    // BUSCO results for MULTIQC
    BUSCO.out.short_summaries_txt
    | ifEmpty ( [ [], [] ] )
    | set { multiqc }


    emit:
    first_table = ch_first_table          // channel: [ val(meta), path(full_table) ] 
    all_tables  = ch_indexed_buscos       // channel: [ val(meta), path(full_tables) ]
    blastp_txt  = DIAMOND_BLASTP.out.txt  // channel: [ val(meta), path(txt) ]
    taxon_id    = ch_taxid                // channel: taxon_id
    multiqc                               // channel: [ meta, summary ]
    versions    = ch_versions             // channel: [ versions.yml ]
}
