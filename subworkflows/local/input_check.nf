//
// Check input samplesheet and get aligned read channels
//

include { samplesheetToList         } from 'plugin/nf-schema'


include { UNTAR                     } from '../../modules/nf-core/untar/main'
include { CAT_CAT                   } from '../../modules/nf-core/cat/cat/main'
include { SAMTOOLS_FLAGSTAT         } from '../../modules/nf-core/samtools/flagstat/main'
include { GENERATE_CONFIG           } from '../../modules/local/generate_config'

workflow INPUT_CHECK {
    take:
    samplesheet     // channel: /path/to/samplesheet
    fasta           // channel: [ meta, path(fasta) ]
    taxon           // channel: val(taxon)
    busco_lin       // channel: val([busco_lin])
    lineage_tax_ids // channel: /path/to/lineage_tax_ids
    databases

    main:
    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Decompress databases if needed
    //

    // Check which need to be decompressed
    ch_dbs_for_untar = databases
        .branch {
            untar: it.size() > 1 && it[1] && it[1].toString().endsWith(".tar.gz")
            skip: true
        }

    // Untar the databases
    UNTAR ( ch_dbs_for_untar.untar )
    ch_versions = ch_versions.mix( UNTAR.out.versions )

    // Join and format dbs
    ch_databases = ch_dbs_for_untar.skip
        .mix( UNTAR.out.untar )
        .map { meta, db -> [ meta + [id: db.baseName], db] }
        .map { db_meta, db_path ->
            if (db_meta.type in ["blastp", "blastx"] && db_path.isDirectory()) {
                [db_meta, file(db_path.toString() + "/${db_path.name}", checkIfExists: true)]
            } else {
                [db_meta, db_path]
            }
        }
        .branch { db_meta, db_path ->
            blastn: db_meta.type == "blastn"
            blastp: db_meta.type == "blastp"
            blastx: db_meta.type == "blastx"
            precomputed_busco: db_meta.type == "precomputed_busco"
            busco: db_meta.type == "busco"
            taxdump: db_meta.type == "taxdump"
        }

    //
    // SUBWORKFLOW: Process samplesheet
    //
    if ( params.fetchngs_samplesheet ) {
        Channel
            .fromList(samplesheetToList(samplesheet, "assets/schema_fetchngs_input.json"))
            .map {it[0]}
            .branch { row ->
                paired: row.fastq_2
                    // Reformat for CAT_CAT
                    [[id: row.run_accession, row:row], [row.fastq_1, row.fastq_2]]
                not_paired: true
            }
            .set { reads_pairedness }

        CAT_CAT ( reads_pairedness.paired )
        ch_versions = ch_versions.mix ( CAT_CAT.out.versions.first() )

        CAT_CAT.out.file_out
        | map { meta, file -> meta.row + [fastq_1: file] }
        | mix ( reads_pairedness.not_paired )
        | map { create_data_channels_from_fetchngs(it) }
        | set { read_files }

    } else {
        Channel
            .fromList(samplesheetToList(samplesheet, "assets/schema_input.json"))
            .map { check_data_channel(it) }
            .set { read_files }
    }


    // Extract the read counts
    SAMTOOLS_FLAGSTAT ( read_files.map { meta, datafile -> [meta, datafile, []] } )
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions.first())

    read_files
    | join( SAMTOOLS_FLAGSTAT.out.flagstat )
    | map { meta, datafile, stats -> [meta + get_read_counts(stats), datafile] }
    | set { reads }


    // Get the source paths of all the databases, except Busco which is not recorded in the blobDir meta.json
    databases
    | filter { meta, file -> meta.type != "busco" && meta.type != "precomputed_busco" }
    | map {meta, file -> [meta, file.toUriString()]}
    | set { db_paths }


    GENERATE_CONFIG (
        fasta,
        taxon,
        busco_lin,
        ch_databases.blastn,
        lineage_tax_ids,
        reads.collect(flat: false).ifEmpty([]),
        db_paths.collect(flat: false),
    )
    ch_versions = ch_versions.mix(GENERATE_CONFIG.out.versions.first())


    //
    // Parse the CSV file
    //
    GENERATE_CONFIG.out.csv
    | map { meta, csv -> csv }
    | splitCsv(header: ['key', 'value'])
    | branch {
        taxon_id: it.key == "taxon_id"
                    return it.value
        busco_lineage: it.key == "busco_lineage"
                    return it.value
    }
    | set { ch_parsed_csv }


    //
    // Get the taxon ID if we do taxon filtering in blast* searches
    //
    ch_parsed_csv.taxon_id
    | map { params.skip_taxon_filtering ? '' : it }
    | first
    | set { ch_taxon_id }


    //
    // Get the BUSCO linages
    //
    ch_parsed_csv.busco_lineage
    | collect
    | set { ch_busco_lineages }

    // Format pre-computed BUSCOs (if provided)
    // Parse the BUSCO output directories
    ch_parsed_busco = ch_databases.precomputed_busco
        .flatMap { meta, dir ->
            def subdirs = file(dir).listFiles().findAll { it.isDirectory() }
            subdirs.collect { subdir ->
                def lineage = subdir.name.split('_')[1..-1].join('_')
                [[type: 'precomputed_busco', id: subdir.name, lineage: lineage], subdir]
            }
        }

    // Remove any invalid lineages from precomputed_busco
    ch_busco_lineages_list = ch_busco_lineages.flatten()
    ch_parsed_busco_filtered = ch_parsed_busco
        .filter { meta, path ->
            ch_busco_lineages.contains(meta.lineage)
        }
    ch_parsed_busco_filtered = ch_parsed_busco_filtered.ifEmpty { Channel.value([]) }

    //
    // Get the BUSCO path if set
    //
    ch_databases.busco
    | map { _, db_path -> db_path }
    | ifEmpty( [] )
    | first
    | set { ch_busco_db }

    emit:
    reads                                   // channel: [ val(meta), path(datafile) ]
    config = GENERATE_CONFIG.out.yaml       // channel: [ val(meta), path(yaml) ]
    synonyms_tsv = GENERATE_CONFIG.out.synonyms_tsv     // channel: [ val(meta), path(tsv) ]
    categories_tsv = GENERATE_CONFIG.out.categories_tsv // channel: [ val(meta), path(tsv) ]
    taxon_id = ch_taxon_id                  // channel: val(taxon_id)
    busco_lineages = ch_busco_lineages      // channel: val([busco_lin])
    blastn = ch_databases.blastn.first()    // channel: [ val(meta), path(blastn_db) ]
    blastp = ch_databases.blastp.first()    // channel: [ val(meta), path(blastp_db) ]
    blastx = ch_databases.blastx.first()    // channel: [ val(meta), path(blastx_db) ]
    precomputed_busco = ch_parsed_busco    // channel: [ val(meta), path(busco_run_dir) ]
    busco_db = ch_busco_db                  // channel: [ path(busco_db) ]
    taxdump = ch_databases.taxdump.map { _, db_path -> db_path }          // channel: [ path(taxdump) ]
    versions = ch_versions                  // channel: [ versions.yml ]
}

// Function to get list of [ meta, datafile ]
def check_data_channel(meta, datafile) {

    if ( !params.align && !datafile.name.endsWith(".bam") && !datafile.name.endsWith(".cram") ) {
        error("ERROR: Please check input samplesheet and pipeline parameters -> Data file is in FastA/FastQ format but --align is not set!\n${datafile}")
    }

    return [meta, datafile]
}

// Function to get list of [ meta, datafile ]
def create_data_channels_from_fetchngs(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.run_accession
    meta.layout     = row.library_layout

    // Same as https://github.com/blobtoolkit/blobtoolkit/blob/4.3.3/src/blobtoolkit-pipeline/src/lib/functions.py#L30-L39
    // with the addition of "hic"
    switch (row.instrument_platform) {
        case "ILLUMINA":
            meta.datatype = (row.library_strategy == "Hi-C" ? "hic" : "illumina")
            break
        case "OXFORD_NANOPORE":
            meta.datatype = "ont"
            break
        case "PACBIO_SMRT":
            meta.datatype = (row.instrument_model == "Sequel" ? "pacbio_clr" : "pacbio")
            break
        default:
            meta.datatype = "illumina"
    }

    return [ meta, file(row.fastq_1) ]
}

// Function to get the read counts from a samtools flagstat file
def get_read_counts ( stats ) {
    // create meta map
    def read_count_meta = [:]

    // Read the first line of the flagstat file
    // 3127898040 + 0 in total (QC-passed reads + QC-failed reads)
    // and make the sum of both integers
    stats.withReader {
        line = it.readLine()
        def lspl = line.split()
        def read_count = lspl[0].toLong() + lspl[2].toLong()
        read_count_meta.read_count = read_count
    }

    return read_count_meta
}
