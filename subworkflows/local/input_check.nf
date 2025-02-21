//
// Check input samplesheet and get aligned read channels
//

include { UNTAR                     } from '../../modules/nf-core/untar/main'
include { CAT_CAT                   } from '../../modules/nf-core/cat/cat/main'
include { SAMTOOLS_FLAGSTAT         } from '../../modules/nf-core/samtools/flagstat/main'
include { SAMPLESHEET_CHECK         } from '../../modules/local/samplesheet_check'
include { FETCHNGSSAMPLESHEET_CHECK } from '../../modules/local/fetchngssamplesheet_check'
include { GENERATE_CONFIG           } from '../../modules/local/generate_config'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    fasta       // channel: [ meta, path(fasta) ]
    taxon       // channel: val(taxon)
    busco_lin   // channel: val([busco_lin])
    lineage_tax_ids        // channel: /path/to/lineage_tax_ids
    blastn                  // channel: [ path(blastn_db) ]
    blastp                  // channel: [ path(blastp_db) ]
    blastx                  // channel: [ path(blastx_db) ]
    busco_output            // channel: [ path(busco_output) ]
    busco_db                // channel: [ path(busco_db) ]
    taxdump                 // channel: [ path(taxdump) ]

    main:
    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Decompress databases if needed
    //

    // Join into single databases channel, filtering out null values
    databases = Channel.empty()
    [blastn, blastp, blastx, busco_db, busco_output, taxdump].each { channel ->
        if (channel) {
            databases = databases.mix(channel.map { it -> it instanceof List ? it : [it] })
        }
    }

    // Check which need to be decompressed, handling potential null values
    ch_dbs_for_untar = databases
        .branch {
            untar: it.size() > 1 && it[1] && it[1].toString().endsWith(".tar.gz")
            skip: true
        }

    // Untar the databases
    UNTAR ( ch_dbs_for_untar.untar.ifEmpty([]) )
    ch_versions = ch_versions.mix( UNTAR.out.versions.ifEmpty([]) )

    // Join and format dbs
    ch_databases = Channel.empty()
    if (UNTAR.out.untar) ch_databases = ch_databases.mix(UNTAR.out.untar)
    if (ch_dbs_for_untar.skip) ch_databases = ch_databases.mix(ch_dbs_for_untar.skip)

    ch_databases = ch_databases
        .filter { it[1] != null }  // Filter out any remaining null values
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
            busco_output: db_meta.type == "busco_output"
            busco: db_meta.type == "busco"
            taxdump: db_meta.type == "taxdump"
        }

    //
    // SUBWORKFLOW: Process samplesheet
    //
    if ( params.fetchngs_samplesheet ) {
        FETCHNGSSAMPLESHEET_CHECK ( samplesheet )
            .csv
            .splitCsv ( header:true, sep:',', quote:'"' )
            .branch { row ->
                paired: row.fastq_2
                    [[id: row.run_accession, row:row], [row.fastq_1, row.fastq_2]]
                not_paired: true
            }
            .set { reads_pairedness }
        ch_versions = ch_versions.mix ( FETCHNGSSAMPLESHEET_CHECK.out.versions.first() )

        CAT_CAT ( reads_pairedness.paired )
        ch_versions = ch_versions.mix ( CAT_CAT.out.versions.first() )

        CAT_CAT.out.file_out
        | map { meta, file -> meta.row + [fastq_1: file] }
        | mix ( reads_pairedness.not_paired )
        | map { create_data_channels_from_fetchngs(it) }
        | set { read_files }

    } else {
        SAMPLESHEET_CHECK ( samplesheet )
            .csv
            .splitCsv ( header:true, sep:',' )
            .map { create_data_channels(it) }
            .set { read_files }
        ch_versions = ch_versions.mix ( SAMPLESHEET_CHECK.out.versions.first() )
    }


    // Extract the read counts
    SAMTOOLS_FLAGSTAT ( read_files.map { meta, datafile -> [meta, datafile, []] } )
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions.first())

    read_files
    | join( SAMTOOLS_FLAGSTAT.out.flagstat )
    | map { meta, datafile, stats -> [meta + get_read_counts(stats), datafile] }
    | set { reads }


    // Format pre-computed BUSCOs (if provided)
    // Parse the BUSCO output directories
    if (ch_databases.busco_output) {
        ch_parsed_busco = ch_databases.busco_output
            .flatMap { meta, dir ->
                def subdirs = file(dir).listFiles().findAll { it.isDirectory() }
                subdirs.collect { subdir ->
                    def lineage = subdir.name.split('_')[1..-1].join('_')
                    [[type: 'busco_output', id: subdir.name, lineage: lineage], subdir]
                }
            }
    } else {
        ch_parsed_busco = Channel.empty()
    }


    GENERATE_CONFIG (
        fasta,
        taxon,
        busco_lin,
        lineage_tax_ids,
        reads.collect(flat: false).ifEmpty([]),
        ch_databases.blastp,
        ch_databases.blastx,
        ch_databases.blastn,
        ch_databases.taxdump,
        ch_parsed_busco.toList()
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

    // Remove any invalid lineages from busco_outputs
    ch_busco_lineages_list = ch_busco_lineages.flatten()
    ch_parsed_busco_filtered = ch_parsed_busco
        .filter { meta, path ->
            ch_busco_lineages.contains(meta.lineage)
        }
    ch_parsed_busco_filtered = ch_parsed_busco_filtered.ifEmpty { Channel.empty() }

    emit:
    reads                                                        // channel: [ val(meta), path(datafile) ]
    config = GENERATE_CONFIG.out.yaml                            // channel: [ val(meta), path(yaml) ]
    synonyms_tsv = GENERATE_CONFIG.out.synonyms_tsv              // channel: [ val(meta), path(tsv) ]
    categories_tsv = GENERATE_CONFIG.out.categories_tsv          // channel: [ val(meta), path(tsv) ]
    taxon_id = ch_taxon_id                                       // channel: val(taxon_id)
    busco_lineages = ch_busco_lineages                           // channel: val([busco_lin])
    blastn = ch_databases.blastn                                 // channel: [ val(meta), path(blastn_db) ]
    blastp = ch_databases.blastp                                 // channel: [ val(meta), path(blastp_db) ]
    blastx = ch_databases.blastx                                 // channel: [ val(meta), path(blastx_db) ]
    busco_output = ch_parsed_busco_filtered                      // channel: [ val(meta), path(busco_output) ]
    busco_db = ch_databases.busco.map { _, db_path -> db_path }  // channel: [ path(busco_db) ]
    taxdump = ch_databases.taxdump.map { _, db_path -> db_path } // channel: [ path(taxdump) ]
    versions = ch_versions                                       // channel: [ versions.yml ]
}

// Function to get list of [ meta, datafile ]
def create_data_channels(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.datatype   = row.datatype
    meta.layout     = row.library_layout

    // add path(s) of the read file(s) to the meta map
    def data_meta = []

    if ( !params.align && !row.datafile.endsWith(".bam") && !row.datafile.endsWith(".cram") ) {
        exit 1, "ERROR: Please check input samplesheet and pipeline parameters -> Data file is in FastA/FastQ format but --align is not set!\n${row.datafile}"
    }

    if ( !file(row.datafile).exists() ) {
        exit 1, "ERROR: Please check input samplesheet -> Data file does not exist!\n${row.datafile}"
    } else {
        data_meta = [ meta, file(row.datafile) ]
    }

    return data_meta
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


    // add path(s) of the read file(s) to the meta map
    def data_meta = []

    if ( !file(row.fastq_1).exists() ) {
        exit 1, "ERROR: Please check input samplesheet -> Data file does not exist!\n${row.fastq_1}"
    } else {
        data_meta = [ meta, file(row.fastq_1) ]
    }

    return data_meta
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
