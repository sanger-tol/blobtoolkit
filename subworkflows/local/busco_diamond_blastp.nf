//
// Run BUSCO for a genome and runs diamond_blastp
//

include { BUSCO_BUSCO               } from '../../modules/nf-core/busco/busco/main'
include { BLOBTOOLKIT_EXTRACTBUSCOS } from '../../modules/local/blobtoolkit/extractbuscos'
include { DIAMOND_BLASTP            } from '../../modules/nf-core/diamond/blastp/main'
include { RESTRUCTUREBUSCODIR       } from '../../modules/sanger-tol/restructurebuscodir/main'


// MODIFY SO ONLY THE OUTPUT OF SANGER-TOL BUSCO CAN BE USED AS INPUT TO THE PIPELINE
// ADD WARNING TO NOTIFY USER

workflow BUSCO_DIAMOND {
    take:
    fasta        // channel: [ val(meta), path(fasta) ]
    busco_lin    // channel: val([busco_lineages])
    busco_db     // channel: path(busco_db)
    odb_version  // channel: val(odb_version)
    blastp       // channel: path(blastp_db)
    taxon_id     // channel: val(taxon_id)
    precomputed_busco // channel: [ val(meta}, path(busco_run_dir) ] optional precomputed busco outputs

    main:
    //
    // LOGIC: Prepare the BUSCO lineages
    //

    // 0. Initialise the basal lineages according to the odb version
    def basal_lineages = [ "eukaryota", "bacteria", "archaea" ]

    ch_basal_lineages = channel.from(basal_lineages)
        .combine(odb_version)
        .map { lineage, version -> lineage + version }

    // Combine the list of relevant lineages with the basal lineages, and with the fasta
    // 1. Convert the list of strings to a channel of a strings
    ch_fasta_with_lineage = busco_lin
        .flatMap()
        // 2. Add the basal lineages, and remove any duplicate introduced
        .concat(ch_basal_lineages)
        .unique()
        // 3. Add a (0-based) index to record the original order (i.e. by age) – withIndex doesn't work on channels
        .toList()
        .flatMap { lineages -> lineages.withIndex() }
        // 4. Add the genome fasta and meta, and keys for the lineage so that we can distinguish the BUSCO jobs and group their outputs later
        .combine(fasta)
        .map { lineage_name, lineage_index, meta, genome -> [meta + [lineage_name: lineage_name, lineage_index: lineage_index], genome] }


    //
    // LOGIC: Format pre-computed outputs
    //
    ch_precomputed_busco = precomputed_busco
        .map { meta, dir -> [meta.lineage, [meta, dir]] }

    ch_combined = ch_fasta_with_lineage
        .map {
            meta, file -> [meta.lineage_name, [meta, file]]
        }
        .join(ch_precomputed_busco, by: 0, remainder: true)
        .map { _lineage, fasta_data, busco_data ->
            def (meta, file) = fasta_data
            def (_busco_meta, busco_dir) = busco_data ?: [null, null]
            [meta + [busco_dir: busco_dir], file]
        }

    // NOTE: Branch based on whether there's a pre-computed BUSCO output
    ch_busco_to_run = ch_combined.branch { meta, _fasta ->
        precomputed: meta.busco_dir != null
        to_compute: true
    }


    //
    // MODULE: Run BUSCO search
    //
    BUSCO_BUSCO(
        ch_busco_to_run.to_compute,
        'genome',
        ch_busco_to_run.to_compute.map { meta, _fasta -> meta.lineage_name },
        busco_db,
        [],
        []
    )


    //
    // LOGIC: Join new and pre-computed BUSCO outputs
    //
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


    //
    // MODULE: Tidy up the BUSCO output directories before publication
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
                    outputs.seq_dir ?: [],
                ]
            }
    )

    //
    // LOGIC: Format pre-computed BUSCO outputs
    //
    ch_formatted_precomputed = ch_busco_to_run.precomputed
        .map { meta, _fasta ->
            def busco_dir = file(meta.busco_dir)
            [
                meta,
                file("${busco_dir}/*.${meta.lineage_name}.multi_copy_busco_sequences.tar.gz"),
                file("${busco_dir}/*.${meta.lineage_name}.single_copy_busco_sequences.tar.gz"),
                file("${busco_dir}/*.${meta.lineage_name}.fragmented_busco_sequences.tar.gz"),
                file("${busco_dir}/*.${meta.lineage_name}.full_table.tsv.gz"),
            ]
        }
        .multiMap{ meta, multi, single, frag, full_table_tsv ->
            busco_seq: tuple(meta, multi, single, frag)
            full_table: tuple(meta, [full_table: full_table_tsv])
        }

    //
    // MODULE: COMBINE THE BUSCO FOLDERS INTO THE EXPECTED FORMAT
    //
    // NOTE: THIS IS SOMEWHAT INEFFICIENT BECAUSE IT IS HAPPENING BEFORE FILTERING
    //       FOR ONLY THE BASAL STUFF SO MOST OF THE OUTPUT IS UNUSED
    ch_precomputed_buscos = PREPARE_PRECOMPUTED_BUSCOS(ch_formatted_precomputed.busco_seq)

    formatted_precomputed_buscos = ch_precomputed_buscos
        .map { meta, targz ->
            [
                meta.lineage_name,
                meta,
                targz,
            ]
        }
        .join(ch_basal_lineages)
        .collect(flat: false)  { _lineage_name, _meta, outputs -> outputs }


    //
    // LOGIC: Select input for BLOBTOOLKIT_EXTRACTBUSCOS
    //
    ch_basal_buscos = ch_all_busco_outputs
        .map { meta, outputs -> [meta.lineage_name, meta, outputs] }
        // The join is equivalent to selecting the channel items whose lineage is basal
        .join(ch_basal_lineages)
        // Without flat:false, collect will flatten meta and outputs
        .collect(flat: false) { _lineage_name, _meta, outputs -> outputs.seq_dir }

    //
    // MODULE: Extract BUSCO genes from the basal lineages
    //
    btk_extract_input = ch_basal_buscos.combine(formatted_precomputed_buscos)
    btk_extract_input.view{"All BUSCO: $it"}
    BLOBTOOLKIT_EXTRACTBUSCOS (
        fasta,
        btk_extract_input
    )


    //
    // LOGIC: Align BUSCO genes against the BLASTp database
    //       FILTER OUT EMPTY FILES FROM ANALYSIS
    //
    ch_busco_genes = BLOBTOOLKIT_EXTRACTBUSCOS.out.genes
        .filter { _meta, file -> file.size() > 140 &&
                    params.blast_annotations in ["all", "blastp", "blastx"]
        }


    //
    // MODULE: Hardcoded to match the format expected by blobtools
    //         DIAMOND WILL NOT RUN IF blast_annotations IS SET TO `off`
    //
    // NOTE:   BLASTP has an issue with sometimes not being deterministic
    //         Matthieu is checking, we can chack as part of this too.
    //
    def outfmt = 6
    def cols   = 'qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
    DIAMOND_BLASTP (
        ch_busco_genes,
        blastp,
        outfmt,
        cols,
        taxon_id
    )


    //
    // MODULE: Order BUSCO results according to the lineage index
    //
    // NOTE: BRANCH FULL TABLES GO HEREE, should be just a mix
    precomputed_busco = ch_formatted_precomputed.full_table

    ch_indexed_buscos = ch_all_busco_outputs
        // 0. Filter out the BUSCO results that found no gene (seen for archaea/bacteria)
        .filter { _meta, outputs -> outputs.full_table }
        // 1. Extract the necessary information and create a consistent structure
        .map { meta, outputs ->
            def cleaned_meta = meta.findAll { pair -> pair.key != "lineage_name" && pair.key != "lineage_index" && pair.key != "busco_dir" }
            def full_table = outputs.full_table
            def lineage_index = meta.lineage_index
            [cleaned_meta, [full_table, lineage_index]]
        }
        // 2. Group by the cleaned meta information
        .groupTuple(by: 0)
        // 3. Sort the tables by lineage index and collect only the tables
        .map { meta, table_positions ->
            [
                meta,
                table_positions.sort { a, b -> a[1] <=> b[1] }.collect { name, _pos -> name }
            ]
        }


    // Select BUSCO results for taxonomically closest database
    ch_first_table = ch_indexed_buscos
        .map { meta, tables -> [meta, tables[0]] }

    // BUSCO results for MULTIQC
    multiqc = BUSCO_BUSCO.out.short_summaries_txt
        .map { _meta, outputs -> outputs }

    emit:
    first_table = ch_first_table          // channel: [ val(meta), path(full_table) ]
    all_tables  = ch_indexed_buscos       // channel: [ val(meta), path(full_tables) ]
    blastp_txt  = DIAMOND_BLASTP.out.txt  // channel: [ val(meta), path(txt) ]
    multiqc                               // channel: [ meta, summary ]
}

process PREPARE_PRECOMPUTED_BUSCOS {
    tag "${meta.id}"
    label 'process_single'
    container "docker.io/genomehubs/blobtoolkit:4.4.6"

    input:
    tuple val(meta), path(multi_copy), path(single_copy), path(fragmented)

    output:
    tuple val(meta), path("busco_sequences.tar.gz"), emit: busco

    when: task.ext.when == null || task.ext.when

    script:
    """
    python - <<'PY'
    import tarfile
    from pathlib import Path

    output = Path("busco_sequences.tar.gz")

    archives = [
        Path("${multi_copy}"),
        Path("${single_copy}"),
        Path("${fragmented}"),
    ]

    # Write to a gz file as destination
    with tarfile.open(output, "w:gz") as dst:
        # For archive INSIDE the larger archive
        for archive in archives:
            # read the gzip file
            with tarfile.open(archive, "r:gz") as src:
                for member in src:
                    # modify the name (member section, e.g. the subfolder) - (group inside extract-busco-genes)
                    # This is essential thanks to some regex inside of the tool
                    member.name = f"busco_sequences/{member.name}"
                    if member.isfile():
                        with src.extractfile(member) as fh:
                            dst.addfile(member, fh)
                    else:
                        dst.addfile(member)
    PY """
}
