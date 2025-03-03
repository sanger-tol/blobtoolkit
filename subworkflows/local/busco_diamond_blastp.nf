//
// Run BUSCO for a genome and runs diamond_blastp
//

include { BUSCO_BUSCO as BUSCO      } from '../../modules/nf-core/busco/busco/main'
include { BLOBTOOLKIT_EXTRACTBUSCOS } from '../../modules/local/blobtoolkit/extractbuscos'
include { DIAMOND_BLASTP            } from '../../modules/nf-core/diamond/blastp/main'
include { RESTRUCTUREBUSCODIR       } from '../../modules/local/restructurebuscodir'


workflow BUSCO_DIAMOND {
    take:
    fasta        // channel: [ val(meta), path(fasta) ]
    busco_lin    // channel: val([busco_lineages])
    busco_db     // channel: path(busco_db)
    blastp       // channel: path(blastp_db)
    taxon_id     // channel: val(taxon_id)


    main:
    ch_versions = Channel.empty()


    //
    // Prepare the BUSCO linages
    //
    // 0. Initialise sone variables
    basal_lineages = [ "eukaryota_odb10", "bacteria_odb10", "archaea_odb10" ]
    def lineage_position = 0
    // 1. Start from the taxon's lineages
    busco_lin
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
    DIAMOND_BLASTP ( ch_busco_genes, blastp, outext, cols, taxon_id )
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
    multiqc                               // channel: [ meta, summary ]
    versions    = ch_versions             // channel: [ versions.yml ]
}
