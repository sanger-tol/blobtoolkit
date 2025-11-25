//
// Check input samplesheet and get aligned read channels
//

include { samplesheetToList         } from 'plugin/nf-schema'
include { UNTAR                     } from '../../modules/nf-core/untar/main'
include { CAT_CAT                   } from '../../modules/nf-core/cat/cat/main'
include { SAMTOOLS_FLAGSTAT         } from '../../modules/nf-core/samtools/flagstat/main'
include { GENERATE_CONFIG           } from '../../modules/local/generate_config'
include { JSONIFY_TAXDUMP           } from '../../modules/local/jsonify_taxdump'

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
    // MODULE: Check which need to be decompressed & Untar if needed
    //
    ch_dbs_for_untar = databases
        .branch { db_meta, db_path ->
            untar: db_path.name.endsWith( ".tar.gz" )
            skip: true
        }

    UNTAR ( ch_dbs_for_untar.untar )
    ch_versions = ch_versions.mix( UNTAR.out.versions.first() )


    //
    // MODULE: Join and format dbs
    //
    ch_databases = ch_dbs_for_untar.skip
        .mix( UNTAR.out.untar )
        .map { meta, db -> [ meta + [id: db.baseName], db] }
        .map { db_meta, db_path ->
            if (db_meta.type in ["blastp", "blastx"] && db_path.isDirectory()) {
                [db_meta, file(db_path.toString() + "/${db_path.name}", checkIfExists: true)]
            } else if (db_meta.type == "blastn") {
                // Special handling for BLAST nucleotide databases
                // If db_path is a directory (from untar), look for .nal file inside
                def actual_db_path = db_path
                if (db_path.isDirectory()) {
                    def nal_files = db_path.listFiles().findAll { it.name.endsWith('.nal') }
                    if (nal_files.size() == 1) {
                        actual_db_path = nal_files[0]
                    } else if (nal_files.size() > 1) {
                        error """
                        ERROR: Multiple .nal files found in extracted blastn database directory: ${db_path}
                        Found: ${nal_files.collect { it.name }.join(', ')}
                        Please ensure the database archive contains only one .nal file.
                        """
                    } else {
                        error """
                        ERROR: No .nal file found in extracted blastn database directory: ${db_path}
                        A valid BLAST nucleotide database must contain a .nal file.
                        """
                    }
                }
                def (resolved_path, db_name) = validateBlastnDatabase(actual_db_path)
                [db_meta, resolved_path]
            } else if (db_meta.type == "busco") {
                // Special handling for BUSCO databases
                def resolved_path = validateBuscoDatabase(db_path)
                [db_meta, resolved_path]
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
    // MODULE: Process samplesheet
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


    //
    // MODULE: Extract the read counts and insert into the meta
    //
    SAMTOOLS_FLAGSTAT (
        read_files.map { meta, datafile -> [meta, datafile, []] }
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions.first())

    read_files
    | join( SAMTOOLS_FLAGSTAT.out.flagstat )
    | map { meta, datafile, stats -> [meta + get_read_counts(stats), datafile] }
    | set { reads }


    //
    // MODULE:  Get the source paths of all the databases, except Busco which is not recorded in the blobDir meta.json
    //
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
    // LOGIC: Parse the CSV file
    //
    GENERATE_CONFIG.out.csv
    | map { meta, csv -> csv }
    | splitCsv(header: ['key', 'value'])
    | branch {
        taxon_id: it.key == "taxon_id"
                    return it.value
        odb_version: it.key == "odb_version"
                    return it.value
        busco_lineage: it.key == "busco_lineage"
                    return it.value
    }
    | set { ch_parsed_csv }


    //
    // LOGIC: Get the taxon ID if we do taxon filtering in blast* searches
    //
    ch_parsed_csv.taxon_id
    | map { params.skip_taxon_filtering ? '' : it }
    | first
    | set { ch_taxon_id }


    //
    // LOGIC: Get the ODB version to use
    //
    ch_parsed_csv.odb_version
    | collect
    | set { ch_odb_version }


    //
    // LOGIC: Get the BUSCO linages
    //
    ch_parsed_csv.busco_lineage
    | collect
    | set { ch_busco_lineages }

    //
    // LOGIC: Format pre-computed BUSCOs (if provided)
    //          Parse the BUSCO output directories
    ch_parsed_busco = ch_databases.precomputed_busco
        .flatMap { meta, dir ->
            def subdirs = file(dir).listFiles().findAll { it.isDirectory() }
            subdirs.collect { subdir ->
                def lineage = subdir.name.startsWith('run_') ? subdir.name.substring(4) : subdir.name
                [[type: 'precomputed_busco', id: subdir.name, lineage: lineage], subdir]
            }
        }

    //
    // LOGIC: Remove any invalid lineages from precomputed_busco
    //
    //ch_busco_lineages_list = ch_busco_lineages.flatten()
    // ch_parsed_busco_filtered = ch_parsed_busco
    //     .filter { meta, path ->
    //         ch_busco_lineages.contains(meta.lineage)
    //     }
    // ch_parsed_busco_filtered = ch_parsed_busco_filtered.ifEmpty { Channel.value([]) }


    //
    // LOGIC: Get the BUSCO path if set
    //
    ch_databases.busco
    | map { _, db_path -> db_path }
    | ifEmpty( [] )
    | set { ch_busco_db }


    //
    // MODULE: Convert the taxdump to a JSON file if there isn't one yet
    //
    ch_databases.taxdump
    | filter { meta, db_path -> ! db_path.isFile() }
    | map { meta, db_path -> [meta, db_path, db_path.listFiles().find { it.getName().endsWith('.json') }] }
    | branch { meta, db_path, json_path ->
        json: json_path
                return [meta, json_path]
        dir: true
                return [meta, db_path]
    }
    | set { taxdump_dirs }


    JSONIFY_TAXDUMP( taxdump_dirs.dir )
    ch_versions = ch_versions.mix(JSONIFY_TAXDUMP.out.versions.first())


    ch_databases.taxdump
    | filter { meta, db_path -> db_path.isFile() }
    | mix ( taxdump_dirs.json )
    | mix( JSONIFY_TAXDUMP.out.json )
    | map { _, db_path -> db_path }
    | set { ch_taxdump }


    emit:
    reads                                   // channel: [ val(meta), path(datafile) ]
    config = GENERATE_CONFIG.out.yaml       // channel: [ val(meta), path(yaml) ]
    synonyms_tsv = GENERATE_CONFIG.out.synonyms_tsv     // channel: [ val(meta), path(tsv) ]
    categories_tsv = GENERATE_CONFIG.out.categories_tsv // channel: [ val(meta), path(tsv) ]
    taxon_id = ch_taxon_id                  // channel: val(taxon_id)
    busco_lineages = ch_busco_lineages      // channel: val([busco_lin])
    odb_version = ch_odb_version            // channel: val(odb_version)
    blastn = ch_databases.blastn.first()    // channel: [ val(meta), path(blastn_db) ]
    blastp = ch_databases.blastp.first()    // channel: [ val(meta), path(blastp_db) ]
    blastx = ch_databases.blastx.first()    // channel: [ val(meta), path(blastx_db) ]
    precomputed_busco = ch_parsed_busco     // channel: [ val(meta), path(busco_run_dir) ]
    busco_db = ch_busco_db.first()          // channel: [ path(busco_db) ]
    taxdump = ch_taxdump.first()            // channel: [ path(taxdump) ]
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

/*
 * Function to validate and resolve BUSCO database paths
 * Handles the common user error of including '/lineages' at the end of the path
 */
def validateBuscoDatabase(db_path) {
    def path_file = file(db_path)
    if (path_file.isDirectory()) {
        // Check if path ends with /lineages and has a parent directory
        if (path_file.name == 'lineages' && path_file.parent != null) {
            def parent_dir = file(path_file.parent)
            log.info "BUSCO path correction: Detected '/lineages' suffix in path"
            log.info "  Original path: ${path_file}"
            log.info "  Corrected path: ${parent_dir}"
            log.info "This prevents the common error where BUSCO tries to use '${path_file}/lineages/lineage_name' instead of '${parent_dir}/lineages/lineage_name'"
            return parent_dir
        }
        // Check if path points to a specific lineage directory (e.g., eukaryota_odb10)
        else if (path_file.name.endsWith('_odb10') && path_file.parent != null) {
            def parent_dir = file(path_file.parent)
            // Check if parent is 'lineages' - if so, we need to go up two levels
            if (parent_dir.name == 'lineages' && parent_dir.parent != null) {
                def busco_root = file(parent_dir.parent)
                log.info "BUSCO path correction: Detected specific lineage directory in path"
                log.info "  Original path: ${path_file} (specific lineage: ${path_file.name})"
                log.info "  Corrected path: ${busco_root}"
                log.info "This prevents the error where BUSCO tries to use a specific lineage directory instead of the root BUSCO database directory"
                log.warn "Use `--busco_lineages ${path_file.name}` to control the lineage"
                return busco_root
            } else {
                error """
                ERROR: Invalid BUSCO lineage directory structure: ${path_file}
                It appears you're pointing to a specific BUSCO lineage directory (${path_file.name}),
                but the expected directory structure is:
                /path/to/busco_downloads/lineages/${path_file.name}/
                Please provide the path to the root BUSCO database directory.
                Example: --busco /path/to/busco_downloads/
                """
            }
        } else {
            // Path looks correct, return as-is
            log.info "Using BUSCO database path: ${path_file}"
            return path_file
        }
    } else {
        error """
        ERROR: Invalid BUSCO database path: ${path_file}
        BUSCO databases must be directories containing the 'lineages/' subdirectory.
        Please ensure the path points to a valid BUSCO database directory.
        Common issues:
        - Path should point to the directory containing 'lineages/' subdirectory
        - Do NOT include '/lineages' at the end of the path
        - Do NOT point to a specific lineage directory (e.g., eukaryota_odb10)
        - BUSCO databases cannot be individual files
        Example: --busco /path/to/busco_downloads/
        NOT: --busco /path/to/busco_downloads/lineages/
        NOT: --busco /path/to/busco_downloads/lineages/eukaryota_odb10/
        """
    }
}

/*
 * Function to validate and resolve BLAST nucleotide database paths
 * Handles both directory paths (for backwards compatibility) and direct .nal file paths
 */
def validateBlastnDatabase(db_path) {
    def path_file = file(db_path)
    // If user provided a file path, require it to be a .nal file
    if (path_file.isFile()) {
        if (!path_file.name.endsWith('.nal')) {
            error """
            ERROR: Invalid BLAST database file: ${path_file}
            The pipeline requires a BLAST nucleotide database file with a .nal extension.
            Please provide the direct path to the .nal file, for example:
                --blastn /path/to/databases/nt.nal
            """
        }

        if (!path_file.exists()) {
            error """
            ERROR: BLAST database file not found: ${path_file}
            Please check that the path is correct and the file exists.
            """
        }

        def parent_dir = file(path_file.parent)
        def db_name = path_file.name.replaceAll('\\.nal$', '')

        // Create an isolated directory inside the pipeline working directory to avoid
        // exposing other databases in the same parent directory.
        def uuid = java.util.UUID.randomUUID().toString()
        def temp_dir = file("${System.getProperty('user.dir')}/.btk_isolated_${db_name}_${uuid}")
        if (!temp_dir.exists()) {
            temp_dir.mkdirs()
        }

        // Ensure the core BLAST database file types exist in the parent directory.
        // NOTE: files do not strictly need to share the same prefix as the .nal file;
        // many installations may name index/sequence/header files independently.
        def required_exts = ['.nin', '.nsq', '.nhr']
        def parent_listing = parent_dir.listFiles()
        def missing = required_exts.findAll { ext ->
            ! parent_listing.any { it.name.endsWith(ext) }
        }
        if (missing && missing.size() > 0) {
            error """
            ERROR: BLAST database appears incomplete in ${parent_dir}
            Missing required file types: ${missing.join(', ')}
            A valid nucleotide BLAST database directory must contain at least one file of each
            of these types (index, sequence, header), for example:
                someprefix.nin
                someprefix.nsq
                someprefix.nhr
            Note: these files do not necessarily have to match the .nal filename prefix.
            Please check the database files and provide the correct .nal path when ready.
            """
        }

        // Collect files to include in the isolated directory. We include any files that
        // either match the database prefix, or are known BLAST DB extensions, plus
        // optional taxonomy helper files.
        def db_files = parent_listing.findAll { f ->
            f.name.startsWith("${db_name}.") ||
            required_exts.any { ext -> f.name.endsWith(ext) } ||
            f.name in ['taxdb.btd', 'taxdb.bti', 'taxonomy4blast.sqlite3'] ||
            f.name.endsWith('.nal')
        }
        db_files.each { source_file ->
            def link_file = file("${temp_dir}/${source_file.name}")
            if (!link_file.exists()) {
                link_file.createLink(source_file)
            }
        }

        log.info "Direct BLAST database file specified: ${path_file}"
        log.info "Database name: ${db_name}"
        log.info "Created isolated directory: ${temp_dir}"
        log.info "This ensures only the specified database is available to BLAST"
        return [temp_dir, db_name]
    }

    // If a directory is provided, instruct the user to pass a .nal file path instead
    else if (path_file.isDirectory()) {
        error """
        ERROR: BLAST database path is a directory: ${path_file}
        This pipeline requires the direct path to the BLAST .nal file.
        Please pass the .nal file to --blastn, for example:
            --blastn /path/to/databases/nt.nal
        """
    }

    // Fallback for non-existent paths or other invalid inputs
    else {
        error """
        ERROR: Invalid BLAST database path: ${path_file}
        Please provide the direct path to a BLAST .nal file, for example:
            --blastn /path/to/databases/nt.nal
        """
    }
}
