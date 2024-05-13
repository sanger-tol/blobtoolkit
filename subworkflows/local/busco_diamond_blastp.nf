//
// Run BUSCO for a genome and runs diamond_blastp
//

include { NCBI_GET_ODB_TAXON        } from '../../modules/local/get_odb_taxon'
include { BUSCO                     } from '../../modules/nf-core/busco/main'
include { BLOBTOOLKIT_EXTRACTBUSCOS } from '../../modules/local/blobtoolkit/extractbuscos'
include { DIAMOND_BLASTP            } from '../../modules/nf-core/diamond/blastp/main'
include { RESTRUCTUREBUSCODIR       } from '../../modules/local/restructurebuscodir'


workflow BUSCO_DIAMOND {
    take:
    fasta        // channel: [ val(meta), path(fasta) ]
    taxon        // channel: val(taxon)
    busco_lin    // channel: val([busco_lin])
    busco_db     // channel: path(busco_db)
    lineage_tax_ids        // channel: /path/to/lineage_tax_ids
    blastp       // channel: path(blastp_db)


    main:
    ch_versions = Channel.empty()


    //
    // Fetch BUSCO lineages for taxon
    //
    NCBI_GET_ODB_TAXON (
        fasta.combine(taxon).map { meta, fasta, taxon -> [ meta, taxon ] },
        busco_lin,
        lineage_tax_ids,
    )
    ch_versions = ch_versions.mix ( NCBI_GET_ODB_TAXON.out.versions.first() )


    //
    // Get NCBI species ID
    //
    NCBI_GET_ODB_TAXON.out.csv
    | map { meta, csv -> csv }
    | splitCsv(header: ['key', 'value'])
    | filter { it.key == "taxon_id" }
    | map { it.value }
    | set { ch_taxid }


    //
    // Prepare the BUSCO linages
    //
    // 0. Initialise sone variables
    basal_lineages = [ "eukaryota_odb10", "bacteria_odb10", "archaea_odb10" ]
    def lineage_position = 0
    // 1. Parse the NCBI_GET_ODB_TAXON output
    NCBI_GET_ODB_TAXON.out.csv
    | map { meta, csv -> csv }
    | splitCsv(header: ['key', 'value'])
    | filter { it.key == "busco_lineage" }
    | map { it.value }
    | collect
    // 2. Add the (missing) basal lineages
    | map { lineages -> (lineages + basal_lineages).unique() }
    | flatten ()
    // 3. Add a (0-based) index to record the original order (i.e. by age)
    | map { lineage_name -> [lineage_name, lineage_position++] }
    // 4. Move the lineage information to `meta` to be able to distinguish the BUSCO jobs and group their outputs later
    | combine ( fasta )
    | map { lineage_name, lineage_index, meta, genome -> [meta + [lineage_name: lineage_name, lineage_index: lineage_index], genome] }
    | set { ch_fasta_with_lineage }


    //
    // Run BUSCO search
    //
    BUSCO (
        ch_fasta_with_lineage,
        "genome",
        ch_fasta_with_lineage.map { it[0].lineage_name },
        busco_db,
        [],
    )
    ch_versions = ch_versions.mix ( BUSCO.out.versions.first() )


    //
    // Tidy up the BUSCO output directories before publication
    //
    RESTRUCTUREBUSCODIR(
        BUSCO.out.seq_dir
        | map { meta, seq -> [meta, meta.lineage_name] }
        | join ( BUSCO.out.batch_summary )
        | join ( BUSCO.out.short_summaries_txt, remainder: true )
        | join ( BUSCO.out.short_summaries_json, remainder: true )
        | join ( BUSCO.out.busco_dir )
        | map { meta, lineage, batch_summary, short_summaries_txt, short_summaries_json, busco_dir -> [meta, lineage, batch_summary, short_summaries_txt ?: [], short_summaries_json ?: [], busco_dir] }
    )
    ch_versions = ch_versions.mix ( RESTRUCTUREBUSCODIR.out.versions.first() )


    //
    // Select input for BLOBTOOLKIT_EXTRACTBUSCOS
    //
    BUSCO.out.seq_dir
    | filter { meta, seq -> basal_lineages.contains(meta.lineage_name) }
    | map { meta, seq -> seq }
    | collect
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

    // Hardcoded to match the format expected by blobtools
    def outext = 'txt'
    def cols   = 'qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
    DIAMOND_BLASTP ( ch_busco_genes, blastp, outext, cols, ch_taxid )
    ch_versions = ch_versions.mix ( DIAMOND_BLASTP.out.versions.first() )


    // Order BUSCO results according to the lineage index
    BUSCO.out.full_table
    // 1. Restore the original meta map, and pull the index as an extra tuple element
    | map { meta, table -> [meta.findAll { it.key != "lineage_name" && it.key != "lineage_index" }, [table, meta.lineage_index]] }
    // 2. Turn to a single-element channel that has the (one and only) meta map, and all the pairs (table, lineage index) concatenated as a list
    | groupTuple()
    // 3. Sort the pairs and discard the index
    | map { meta, table_positions -> [ meta, table_positions.sort { a, b -> a[1] <=> b[1] } . collect { table, lineage_index -> table } ] }
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
