//
// Run BUSCO for a genome from GOAT and runs diamond_blastp
//

include { GOAT_TAXONSEARCH          } from '../../modules/nf-core/goat/taxonsearch/main'
include { BUSCO                     } from '../../modules/nf-core/busco/main'
include { TARGZ                     } from '../../modules/local/targz'
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
    // Run BUSCO search
    //
    GOAT_TAXONSEARCH.out.taxonsearch
    | map { meta, csv -> csv.splitCsv(header:true, sep:'\t', strip:true) }
    | map { row -> row.odb10_lineage.findAll { it != "" } }
    | map { lineages -> [ lineages + [ "bacteria_odb10", "archaea_odb10" ] ] }
    | flatten ()
    | set { ch_lineages }

    BUSCO ( fasta, ch_lineages, busco_db.collect().ifEmpty([]), [] )
    ch_versions = ch_versions.mix ( BUSCO.out.versions.first() )


    //
    // Create input for BLOBTOOLKIT_EXTRACTBUSCOS
    //
    TARGZ ( BUSCO.out.seq_dir )
    ch_versions = ch_versions.mix ( TARGZ.out.versions.first() )
    
    BUSCO.out.full_table
    | join ( TARGZ.out.archive )
    | map { meta, table, seq -> [ [ "id": table.parent.baseName ], table, seq ] }
    | branch {
        meta, table, seq ->
            archaea   : meta.id == "run_archaea_odb10"
            bacteria  : meta.id == "run_bacteria_odb10"
            eukaryota : meta.id == "run_eukaryota_odb10"
    }
    | set { ch_busco }


    // Extract BUSCO genes from the 3 kingdoms
    BLOBTOOLKIT_EXTRACTBUSCOS ( fasta, ch_busco.archaea, ch_busco.bacteria, ch_busco.eukaryota  )
    ch_versions = ch_versions.mix ( BLOBTOOLKIT_EXTRACTBUSCOS.out.versions.first() )


    //
    // Align BUSCO genes against the BLASTp database
    //    
    BLOBTOOLKIT_EXTRACTBUSCOS.out.genes
    | filter { it[1].size() > 140 }
    | set { ch_busco_genes }

    DIAMOND_BLASTP ( ch_busco_genes, blastp, outext, cols )
    ch_versions = ch_versions.mix ( DIAMOND_BLASTP.out.versions.first() )


    // Select BUSCO results for taxonomically closest database
    BUSCO.out.full_table
    | combine ( ch_lineages.toList().map { it[0] } )
    | filter { meta, table, lineage -> table =~ /$lineage/ }
    | map { meta, table, lineage -> [ meta, table ] }
    | set { ch_first_table }


    emit:
    first_table = ch_first_table          // channel: [ val(meta), path(full_table) ] 
    full_table  = BUSCO.out.full_table    // channel: [ val(meta), path(full_tables) ]
    blastp_txt  = DIAMOND_BLASTP.out.txt  // channel: [ val(meta), path(txt) ]
    versions    = ch_versions             // channel: [ versions.yml ]
}
