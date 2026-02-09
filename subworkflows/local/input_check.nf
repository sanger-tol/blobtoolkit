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
                // If db_path is a directory from untar, look for a .nal or .nin file to use as input
                def actual_db_path = db_path
                if (db_path.isDirectory()) {
                    def nal_files = findAllInDir(db_path) { it.name.endsWith('.nal') }
                    def nin_files = findAllInDir(db_path) { it.name.endsWith('.nin') }
                    if (nal_files.size() == 1) {
                        actual_db_path = nal_files[0]
                    } else if (nal_files.size() > 1) {
                        error """
                        ERROR: Multiple .nal files found in blastn database directory: ${db_path}
                        Found: ${nal_files.collect { it.name }.join(', ')}
                        Please ensure the directory contains only one .nal file.
                        """
                    } else if (nin_files.size() == 1) {
                        actual_db_path = nin_files[0]
                    } else if (nin_files.size() > 1) {
                        error """
                        ERROR: Multiple .nin files found in blastn database directory: ${db_path}
                        Found: ${nin_files.collect { it.name }.join(', ')}
                        Please ensure the directory contains only one .nin file.
                        """
                    } else {
                        error """
                        ERROR: No .nal or .nin file found in blastn database directory: ${db_path}
                        Please ensure the directory contains a .nal or .nin file.
                        """
                    }
                }
                def (resolved_path, db_name) = validateBlastnDatabase(actual_db_path)
                [db_meta + [db_name: db_name], resolved_path]
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
            def subdirs = findAllInDir(dir) { it.isDirectory() }
            subdirs.collect { subdir ->
                def lineage = subdir.name.startsWith('run_') ? subdir.name.substring(4) : subdir.name
                [[type: 'precomputed_busco', id: subdir.name, lineage: lineage], subdir]
            }
        }


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
    | map { meta, db_path -> [meta, db_path, findInDir(db_path) { it.name.endsWith('.json') }] }
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


// Helpers to find files in a directory matching a closure

// Return all as a list, potentially empty
def findAllInDir(path, closure) {
    // listFiles() doesn't work on a symlink
    def files = path.toRealPath().listFiles()
    if (!files) return []
    return files.findAll(closure)
}

// Return one or null
def findInDir(path, closure) {
    def files = findAllInDir(path, closure)
    return files ? files[0] : null
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
        // Check if path points to a specific lineage directory (e.g., eukaryota_odb10, eukaryota_odb12)
        // Accept any lineage name that ends with _odb<digits> (odb10, odb12, etc.)
        else if (path_file.name ==~ /.*_odb\d+$/ && path_file.parent != null) {
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
            def lineages_subdir = findInDir(path_file) { it.name == "lineages" }
            if (lineages_subdir) {
                // Path looks correct, return as-is
                return path_file
            }
            error """
            ERROR: Invalid BUSCO lineage directory structure: ${path_file}
            It appears you're pointing to a directory that does not contain a 'lineages' sub-directory
            Please provide the path to the root BUSCO database directory.
            Example: --busco /path/to/busco_downloads/
            """
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
        NOT: --busco /path/to/busco_downloads/lineages/eukaryota_odb10/ (or other lineage directories such as eukaryota_odb12/)
        """
    }
}

/*
 * Function to validate and resolve BLAST nucleotide database paths
 */
def validateBlastnDatabase(path_file) {
    def parent_dir = null
    def db_name = null
    def nal_files = []
    def nin_files = []

    if (path_file.isDirectory()) {
        error """
        ERROR: Invalid BLAST database path: ${path_file}
        Please provide a direct path to a .nal/.nin file.
        """
    } else if (path_file.isFile()) {
        parent_dir = file(path_file.parent)
        if (!path_file.name.endsWith('.nal') && !path_file.name.endsWith('.nin')) {
            error """
            ERROR: Invalid BLAST database file: ${path_file}
            Please provide a direct path to a .nal or .nin file.
            """
        }
        db_name = path_file.name.replaceAll('\\.(nal|nin)$', '')
    } else {
        error """
        ERROR: Invalid BLAST database path: ${path_file}
        Please provide direct path to a .nal/.nin file.
        """
    }
    def alldb_files = findAllInDir(parent_dir) { it.isFile() } .collect { it.name }
    //log.info "BLAST DB directory file names (${parent_dir}): ${alldb_files}"

    // The specific extensions/names we are looking for
    def requiredExtensions = [".nin", ".nhr", ".nsq"]
    def requiredFullFiles = ["taxonomy4blast.sqlite3", "taxdb.btd", "taxdb.bti"]


    nal_files = alldb_files.findAll { it.endsWith('.nal') }
    nin_files = alldb_files.findAll { it.endsWith('.nin') }


    def missingExtensions = requiredExtensions.findAll { ext ->
        !alldb_files.any { it.endsWith(ext) }
    }
    def missingFullFiles = requiredFullFiles.findAll { name ->
        !alldb_files.contains(name)
    }

    def requiredPrefixed = [".nin", ".nhr", ".nsq"]
    def missingPrefixed = requiredPrefixed.findAll { ext ->
        !alldb_files.any { name ->
            name.startsWith("${db_name}") && name.endsWith(ext)
        }
    }

    if (missingExtensions || missingFullFiles || missingPrefixed) {
        def missing_parts = []
        missingExtensions.each { missing_parts << "*${it}" }
        missingFullFiles.each { missing_parts << it }
        missingPrefixed.each { missing_parts << "${db_name}*${it}" }
        error """
        ERROR: BLAST database appears incomplete in ${parent_dir}
        Missing required files: ${missing_parts.join(', ')}
        A valid nucleotide BLAST database must contain:
            .nin (index file)
            .nhr (header file)
            .nsq (sequence file)
            taxonomy4blast.sqlite3 (taxonomy information file)
            taxdb.btd (taxonomy index file)
            taxdb.bti (taxonomy index file)
        """
    }

    log.info "BLAST DB name: ${db_name} in ${parent_dir}"
    log.info "--- Database Integrity Verified: Ready to run BLAST ---"
    return [parent_dir, db_name]
}
