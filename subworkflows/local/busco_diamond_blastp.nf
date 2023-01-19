//
//
// Run BUSCO for a genome from GOAT and runs diamond_blastp
//
//

nextflow.enable.dsl = 2

include { GOAT_TAXONSEARCH    } from '../../modules/nf-core/goat/taxonsearch/main'
include { BUSCO               } from '../../modules/nf-core/busco/main'
include { TAR                 } from '../../modules/local/tar'
include { EXTRACT_BUSCO_GENES } from '../../modules/local/extract_busco_genes'
include { DIAMOND_BLASTP      } from '../../modules/nf-core/diamond/blastp/main'
include { FASTAWINDOWS        } from '../../modules/nf-core/fastawindows/main'
include { CREATE_BED          } from '../../modules/local/create_bed'
include { COUNT_BUSCOGENES    } from '../../modules/local/count_buscogenes'

workflow BUSCO_DIAMOND {
    take:

    //  Tuple [meta, fasta]:
    fasta

    main:

    ch_versions = Channel.empty()


    //
    // Fetch BUSCO lineages for taxon (or taxa)
    //

    GOAT_TAXONSEARCH (
    fasta.map { fa -> [fa[0], "${params.taxon}", "${params.taxa_file}" ? file("${params.taxa_file}") : []] }
    )
    ch_versions = ch_versions.mix(GOAT_TAXONSEARCH.out.versions)

    //
    // Run BUSCO search
    //

    // Make a new channel in which each output file is replaced by its lineages column, with meta propagated
    ch_lineages_goat = GOAT_TAXONSEARCH.out.taxonsearch \
        | map { meta, csv -> [ meta, csv.splitCsv(header:true, sep:'\t', strip:true) ] } \
        | map { meta, row -> [ meta, row.odb10_lineage.findAll { it != '' } ] }    
    //  Channel containing BUSCO lineages for bacteria and archaea
    ch_lineages_prok = fasta.map { fa -> [fa[0], ['bacteria_odb10','archaea_odb10']] }
    // Lineages from GOAT+bacteria+archaea: emits tuples with [meta, lineage]
    ch_lineages = ch_lineages_goat.combine( ch_lineages_prok, by:0 ) \
        | map { id,goat,prok -> [id, goat.plus(prok)] } \
        | transpose()
    // Cross-product of both channels, using meta as the key
    ch_busco_inputs = fasta.combine(ch_lineages, by: 0)

    BUSCO (
    ch_busco_inputs.map { [ it[0], it[1]] },
    ch_busco_inputs.map { it[2] },
    "${params.busco_lineages_path}",  // Please pass this option. We don't want to download the lineage data every time.
    [] // No config
    )

    //look into adding busco_lineage path to meta - create new channel as alexandar uses busco
    ch_versions = ch_versions.mix(BUSCO.out.versions)

    ch_busco = BUSCO.out.busco_dir
                | map { meta, dir ->
                    lin = dir.name.tokenize("-")[1]
                    [ meta, lin, dir ] }

    ch_tsv_path = GrabFiles(ch_busco).groupTuple(by: [0]).view()
        
    // Generate BED File
    FASTAWINDOWS(fasta)
    ch_versions = ch_versions.mix(FASTAWINDOWS.out.versions)

    CREATE_BED(FASTAWINDOWS.out.mononuc)
    ch_versions = ch_versions.mix(CREATE_BED.out.versions)
    
    ch_bed = CREATE_BED.out.bed.view()

    COUNT_BUSCOGENES(ch_tsv_path, ch_bed)

    //channel.fromPath("${BUSCO.out.busco_dir}".map{it[1] + "/**/full_table.tsv"}).view()
    
    // /busco_output/**/full_table.tsv").collectFile($"it".name).view()
    //try with {}
    //try putting glob pattern in fromPath
    //use collect with file path 
    //paramsDir = BUSCO.out.busco_dir.map{it[1]}.view()

    // /*/*/full_table.tsv
    //ch_path = Channel.watchPath("${paramsDir}/*/*/*.tsv").view()

    //    | map {meta, p -> [p + "/*/*"]} \
        //| collectFile()
        //| view()
    //Channel.fromPath("$ch_tsv_path/full_table.tsv") \
        //| collectFile {} \
    //    | view 

    //channel.fromPath("BUSCO.out.busco_dir").view() 
    //dir = BUSCO.out.busco_dir.map{meta, p -> [p + "/*/*/full_table.tsv"]}.view()
    //dir.combine(/*/*/*/full_table.tsv)
    //.view()
    //try removing tsv from param and add it to dirpath
    //params.dir = BUSCO.out.busco_dir.map{meta, p -> [p + "/*/*"]}.view()
    //dirpath = Channel.fromPath("${params.dir}/full_table.tsv").collectFile().view()
    //dirpath = channel.fromPath("${params.dir}/full_table.tsv").collectFile().view()
    //dirpath = channel.fromPath("BUSCO.out.busco_dir.map{meta, p -> [p + '/*/*/full_table.tsv']}.view()")
    //dirpath.view()
        //.collectFile
    //dir.view("$it")

    //channel.fromPath("$baseDir").view()

    /* 
    //
    // Extract BUSCO genes
    //

    // gets file name which is a folder name in paths to busco full tables
    fasta_filename = fasta.map { meta,fa -> [meta, fa.name] }

    // channel: emits paths to busco results for each lineage
    dir = BUSCO.out.busco_dir.combine(fasta_filename, by:0)

    // filter paths to busco full tables for archaea, bacteria and eukaryota
    tbl_a = dir.filter { "$it" =~ /archaea_odb10/ }.map { meta,a,f -> [meta, "$a/$f/run_archaea_odb10/full_table.tsv", "$a/$f/run_archaea_odb10/archaea_odb10_full_table.tsv"] }.collect()
    tbl_b = dir.filter { "$it" =~ /bacteria_odb10/ }.map { meta,b,f -> [meta, "$b/$f/run_bacteria_odb10/full_table.tsv", "$b/$f/run_bacteria_odb10/bacteria_odb10_full_table.tsv"] }.collect()
    tbl_e = dir.filter { "$it" =~ /eukaryota_odb10/ }.map { meta,e,f -> [meta, "$e/$f/run_eukaryota_odb10/full_table.tsv", "$e/$f/run_eukaryota_odb10/eukaryota_odb10_full_table.tsv"] }.collect()
    // create copies of full tables with a lineage identifier, avoids file name collision
    tbl_a = tbl_a.map { meta,t,u -> [meta,file("$t").copyTo("$u")] }.collect()
    tbl_b = tbl_b.map { meta,t,u -> [meta,file("$t").copyTo("$u")] }.collect()
    tbl_e = tbl_e.map { meta,t,u -> [meta,file("$t").copyTo("$u")] }.collect()
    // combine all three channels into a single channel: tuple( meta, a, b, e )
    tbl_ab = tbl_a.combine(tbl_b, by:0)
    tbl_abe = tbl_ab.combine(tbl_e, by:0)

    // filter paths to busco_sequences for archaea, bacteria and eukaryota
    seq_a = dir.filter { "$it" =~ /archaea_odb10/ }.map { meta,a,f -> [meta, "$a/$f/run_archaea_odb10/busco_sequences", "$a/$f/run_archaea_odb10/archaea_odb10_busco_sequences"] }.collect()
    seq_b = dir.filter { "$it" =~ /bacteria_odb10/ }.map { meta,b,f -> [meta, "$b/$f/run_bacteria_odb10/busco_sequences", "$b/$f/run_bacteria_odb10/bacteria_odb10_busco_sequences"] }.collect()
    seq_e = dir.filter { "$it" =~ /eukaryota_odb10/ }.map { meta,e,f -> [meta, "$e/$f/run_eukaryota_odb10/busco_sequences", "$e/$f/run_eukaryota_odb10/eukaryota_odb10_busco_sequences"] }.collect()
    // create copies of busco_sequences with a lineage identifier, avoids file name collision
    seq_a = seq_a.map { meta,t,u -> [meta,file("$t").copyTo("$u")] }.collect()
    seq_b = seq_b.map { meta,t,u -> [meta,file("$t").copyTo("$u")] }.collect()
    seq_e = seq_e.map { meta,t,u -> [meta,file("$t").copyTo("$u")] }.collect()
    // combine all three channels into a single channel: tuple( meta, a, b, e )
    seq_ab = seq_a.combine(seq_b, by:0)
    seq_abe = seq_ab.combine(seq_e, by:0)

    // module: creates input paths for EXTRACT_BUSCO_GENES
    TAR (
    tbl_abe,
    seq_abe
    )

    // module: extract busco genes
    EXTRACT_BUSCO_GENES (
    TAR.out.dir_abe
    )
    ch_versions = ch_versions.mix(EXTRACT_BUSCO_GENES.out.versions)

    //
    // Runs diamond_blastp with the extracted busco genes
    //


    // path to diamond db
    blastp_db = Channel.fromPath(params.diamondblastp_db)

    DIAMOND_BLASTP (
    EXTRACT_BUSCO_GENES.out.fasta,
    blastp_db,
    "${params.blastp_outext}",
    "${params.blastp_cols}"
    )
    ch_versions = ch_versions.mix(DIAMOND_BLASTP.out.versions)

    */

    emit:

    // diamond_blastp outputs
    //txt      = DIAMOND_BLASTP.out.txt

    // BUSCO outputs
    summary = BUSCO.out.batch_summary
    busco_dir = BUSCO.out.busco_dir

    // tool versions
    versions = ch_versions
}

process GrabFiles {
    tag "${meta.id}"
    //tag "${meta.lineage}"
    executor 'local'

    input:
    tuple val(meta), val(lineage) ,path("in")

    output:
    tuple val(meta), path("in/**/${lineage}_full_table.tsv")

    script:
    //def org = "${meta.id}-${meta.lineage}-busco/${meta.id}.fasta/run_${meta.lineage}/full_table.tsv"
    //def renamed = "${meta.id}-${meta.lineage}-busco/${meta.id}.fasta/run_${meta.lineage}/${meta.lineage}_full_table.tsv"
    def org = "in/${meta.id}.fasta/run_${lineage}/full_table.tsv"
    def renamed = "in/${meta.id}.fasta/run_${lineage}/${lineage}_full_table.tsv"
    """
    mv ${org} ${renamed}
    """
}