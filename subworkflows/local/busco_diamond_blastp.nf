//
// Run BUSCO for a genome and runs diamond_blastp
//

include { BUSCO_BUSCO               } from '../../modules/nf-core/busco/busco/main'
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
    precomputed_busco // channel: [ val(meta}, path(busco_run_dir) ] optional precomputed busco outputs

    main:
    ch_versions = Channel.empty()


    //
    // Prepare the BUSCO lineages
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
    // Format pre-computed outputs
    //
    if (params.precomputed_busco){
        ch_precomputed_busco = precomputed_busco
            .map { meta, dir -> [meta.lineage, [meta, dir]] }
    }else{
        ch_precomputed_busco = Channel.empty()
    }

    ch_combined = ch_fasta_with_lineage
        .map { meta, fasta -> [meta.lineage_name, [meta, fasta]] }
        .join(ch_precomputed_busco, by: 0, remainder: true)
        .map { lineage, fasta_data, busco_data ->
            def (meta, fasta) = fasta_data
            def (busco_meta, busco_dir) = busco_data ?: [null, null]
            [meta + [busco_dir: busco_dir], fasta]
        }

    // Branch based on whether there's a pre-computed BUSCO output
    ch_busco_to_run = ch_combined.branch {
        precomputed: it[0].busco_dir != null
        to_compute: true
    }


    // Format precomputed BUSCO outputs
    if (ch_busco_to_run.precomputed) {
        ch_formatted_precomputed = ch_busco_to_run.precomputed
            .map { meta, fasta ->
                def busco_dir = file(meta.busco_dir)
                def lineage = meta.lineage_name
                [
                    meta,
                    [
                        batch_summary: [],
                        short_summaries_txt: file("${busco_dir}/short_summary.txt"),
                        short_summaries_json: file("${busco_dir}/short_summary.json"),
                        full_table: file("${busco_dir}/full_table.tsv"),
                        missing_busco_list: file("${busco_dir}/missing_busco_list.tsv"),
                        single_copy_proteins: file("${busco_dir}/single_copy_proteins.faa"),
                        seq_dir: file("${busco_dir}/busco_sequences"),
                        translated_dir: file("${busco_dir}/translated_proteins"),
                        busco_dir: busco_dir
                    ]
                ]
            }

    } else {
        ch_formatted_precomputed = Channel.empty()
    }

    //
    // Run BUSCO search
    //
    BUSCO_BUSCO(
        ch_busco_to_run.to_compute,
        'genome',
        ch_busco_to_run.to_compute.map { it[0].lineage_name },
        busco_db,
        []
    )
    ch_versions = ch_versions.mix ( BUSCO_BUSCO.out.versions.first() )

    // Join new and pre-computed BUSCO outputs
    ch_all_busco_outputs = BUSCO_BUSCO.out.batch_summary
        .join(BUSCO_BUSCO.out.short_summaries_txt, by: 0, remainder: true )
        .join(BUSCO_BUSCO.out.short_summaries_json, by: 0, remainder: true )
        .join(BUSCO_BUSCO.out.full_table, by: 0, remainder: true )
        .join(BUSCO_BUSCO.out.missing_busco_list, by: 0, remainder: true )
        .join(BUSCO_BUSCO.out.single_copy_proteins, by: 0, remainder: true)
        .join(BUSCO_BUSCO.out.seq_dir, by: 0)
        .join(BUSCO_BUSCO.out.translated_dir, by: 0, remainder: true )
        .join(BUSCO_BUSCO.out.busco_dir, by: 0)
        .map { meta, batch_summary, short_summaries_txt, short_summaries_json, full_table, missing_busco_list, single_copy_proteins, seq_dir, translated_dir, busco_dir ->
            [
                meta,
                [
                    batch_summary: batch_summary,
                    short_summaries_txt: short_summaries_txt,
                    short_summaries_json: short_summaries_json,
                    full_table: full_table,
                    missing_busco_list: missing_busco_list,
                    single_copy_proteins: single_copy_proteins,
                    seq_dir: seq_dir,
                    translated_dir: translated_dir,
                    busco_dir: busco_dir
                ]
            ]
        }
        .mix(ch_formatted_precomputed)

    //
    // Tidy up the BUSCO output directories before publication
    //
    RESTRUCTUREBUSCODIR(
        ch_all_busco_outputs
            .map { meta, outputs ->
                [
                    meta,
                    meta.lineage_name,
                    outputs.batch_summary ?: [],
                    outputs.short_summaries_txt ?: [],
                    outputs.short_summaries_json ?: [],
                    outputs.full_table ?: [],
                    outputs.missing_busco_list ?: [],
                    outputs.seq_dir ? "${outputs.seq_dir}/single_copy_busco_sequences" : [],
                    outputs.seq_dir ? "${outputs.seq_dir}/multi_copy_busco_sequences" : [],
                    outputs.seq_dir ? "${outputs.seq_dir}/fragmented_busco_sequences" : [],
                    outputs.busco_dir ? "${outputs.busco_dir}/hmmer_output" : []
                ]
            }
    )

    //
    // Select input for BLOBTOOLKIT_EXTRACTBUSCOS
    //
    ch_all_busco_outputs
        .filter { meta, outputs -> basal_lineages.contains(meta.lineage_name) }
        .map { meta, outputs -> [meta, outputs.seq_dir] }
        .collect { it[1] }
        .set { ch_basal_buscos }

    // Extract BUSCO genes from the basal lineages
    BLOBTOOLKIT_EXTRACTBUSCOS ( fasta, ch_basal_buscos )
    ch_versions = ch_versions.mix ( BLOBTOOLKIT_EXTRACTBUSCOS.out.versions.first() )

    //
    // Align BUSCO genes against the BLASTp database
    //
    BLOBTOOLKIT_EXTRACTBUSCOS.out.genes
        .filter { it[1].size() > 140 }
        .set { ch_busco_genes }

    // Hardcoded to match the format expected by blobtools
    def outext = 'txt'
    def cols   = 'qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
    DIAMOND_BLASTP ( ch_busco_genes, blastp, outext, cols, taxon_id )
    ch_versions = ch_versions.mix ( DIAMOND_BLASTP.out.versions.first() )

    // Order BUSCO results according to the lineage index
    ch_all_busco_outputs
        // 1. Extract the necessary information and create a consistent structure
        | map { meta, outputs ->
            def cleaned_meta = meta.findAll { it.key != "lineage_name" && it.key != "lineage_index" && it.key != "busco_dir" }
            def full_table = outputs.full_table
            def lineage_index = meta.lineage_index
            [cleaned_meta, [full_table, lineage_index]]
        }
        // 2. Group by the cleaned meta information
        | groupTuple(by: 0)
        // 3. Sort the tables by lineage index and collect only the tables
        | map { meta, table_positions ->
            [
                meta,
                table_positions.sort { a, b -> a[1] <=> b[1] }.collect { it[0] }
            ]
        }
        | set { ch_indexed_buscos }

    // Select BUSCO results for taxonomically closest database
    ch_indexed_buscos
        .map { meta, tables -> [meta, tables[0]] }
        .set { ch_first_table }

    // BUSCO results for MULTIQC
    ch_all_busco_outputs
        .map { meta, outputs -> outputs.batch_summary }
        .ifEmpty ( [ [], [] ] )
        .set { multiqc }

    emit:
    first_table = ch_first_table          // channel: [ val(meta), path(full_table) ]
    all_tables  = ch_indexed_buscos       // channel: [ val(meta), path(full_tables) ]
    blastp_txt  = DIAMOND_BLASTP.out.txt  // channel: [ val(meta), path(txt) ]
    multiqc                               // channel: [ meta, summary ]
    versions    = ch_versions             // channel: [ versions.yml ]
}
