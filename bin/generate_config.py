#!/usr/bin/env python3

import argparse
import dataclasses
import os
import sys
import sqlite3
import requests
import typing
import yaml

NCBI_TAXONOMY_API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/%s"
GOAT_LOOKUP_API = "https://goat.genomehubs.org/api/v2/lookup?searchTerm=%s&result=taxon&size=10&taxonomy=ncbi"
GOAT_RECORD_API = "https://goat.genomehubs.org/api/v2/record?recordId=%s&result=taxon&size=10&taxonomy=ncbi"
NCBI_DATASETS_API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/%s/dataset_report?filters.assembly_version=all_assemblies"
NCBI_SEQUENCE_API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/%s/sequence_reports"

RANKS = [
    "genus",
    "family",
    "order",
    "class",
    "phylum",
    "kingdom",
    "superkingdom",
]

BUSCO_BASAL_LINEAGES = ["eukaryota_odb10", "bacteria_odb10", "archaea_odb10"]


def parse_args(args=None):
    Description = "Produce the various configuration files needed within the pipeline"

    parser = argparse.ArgumentParser(description=Description)
    parser.add_argument("--fasta", required=True, help="Path to the Fasta file of the assembly.")
    parser.add_argument("--taxon_query", required=True, help="Query string/integer for this taxon.")
    parser.add_argument("--lineage_tax_ids", required=True, help="Mapping between BUSCO lineages and taxon IDs.")
    parser.add_argument("--output_prefix", required=True, help="Prefix to name the output files.")
    parser.add_argument("--accession", help="Accession number of the assembly (optional).", default=None)
    parser.add_argument("--busco", help="Requested BUSCO lineages.", default=None)
    parser.add_argument("--nt", required=True, help="Path to the NT database")
    parser.add_argument("--read_id", action="append", help="ID of a read set")
    parser.add_argument("--read_type", action="append", help="Type of a read set")
    parser.add_argument("--read_layout", action="append", help="Layout of a read set")
    parser.add_argument("--read_path", action="append", help="Path of a read set")
    parser.add_argument("--blastp", help="Path to the blastp database", required=True)
    parser.add_argument("--blastx", help="Path to the blastx database", required=True)
    parser.add_argument("--blastn", help="Path to the blastn database", required=True)
    parser.add_argument("--taxdump", help="Path to the taxonomy database", required=True)
    parser.add_argument("--version", action="version", version="%(prog)s 2.0")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        os.makedirs(path, exist_ok=True)


@dataclasses.dataclass
class TaxonInfo:
    taxon_id: int
    organism_name: str
    rank: typing.Optional[str]
    lineage: typing.List["TaxonInfo"]


def make_taxon_info_from_goat(body: dict[str, str]) -> TaxonInfo:
    rank = body["taxon_rank"]
    if rank == "no rank":
        rank = None
    if "lineage" in body:
        lineage = [make_taxon_info_from_goat(b) for b in body["lineage"]]
    else:
        lineage = []
    return TaxonInfo(int(body["taxon_id"]), body["scientific_name"], rank, lineage)


def fetch_taxon_info_from_goat(taxon_name: typing.Union[str, int]) -> TaxonInfo:
    if isinstance(taxon_name, int):
        taxon_id = taxon_name
        record_id = "taxon-%d" % taxon_name
    else:
        # Resolve the taxon_id of the species
        response = requests.get(GOAT_LOOKUP_API % taxon_name).json()
        taxon_id = int(response["results"][0]["result"]["taxon_id"])
        record_id = response["results"][0]["id"]

    # Using API, get the taxon_ids of the species and all parents
    response = requests.get(GOAT_RECORD_API % record_id).json()
    body = response["records"][0]["record"]
    return make_taxon_info_from_goat(body)


def fetch_taxon_info_from_ncbi(taxon_name: typing.Union[str, int], with_lineage=True) -> typing.Optional[TaxonInfo]:
    # Using API, get the taxon_ids of the species and all parents
    response = requests.get(NCBI_TAXONOMY_API % taxon_name).json()
    if "taxonomy" in response["taxonomy_nodes"][0]:
        body = response["taxonomy_nodes"][0]["taxonomy"]
        if with_lineage:
            lineage = [
                fetch_taxon_info_from_ncbi(t, with_lineage=False) for t in reversed(body["lineage"][2:])
            ]  # skip root and cellular_organisms
        else:
            lineage = []
        return TaxonInfo(body["tax_id"], body["organism_name"], body.get("rank"), lineage)


def fetch_taxon_info(taxon_name: typing.Union[str, int]) -> TaxonInfo:
    return fetch_taxon_info_from_ncbi(taxon_name) or fetch_taxon_info_from_goat(taxon_name)


def get_classification(taxon_info: TaxonInfo) -> typing.Dict[str, str]:
    ancestors: typing.Dict[str, str] = {}
    for anc_taxon_info in taxon_info.lineage:
        if anc_taxon_info.rank:
            ancestors[anc_taxon_info.rank.lower()] = anc_taxon_info.organism_name
    # https://ncbiinsights.ncbi.nlm.nih.gov/2024/06/04/changes-ncbi-taxonomy-classifications/
    # "superkingdom" will be called "domain"
    if "superkingdom" not in ancestors:
        ancestors["superkingdom"] = ancestors["domain"]
    return {r: ancestors[r] for r in RANKS if r in ancestors}


def get_odb(taxon_info: TaxonInfo, lineage_tax_ids: str, requested_buscos: typing.Optional[str]) -> typing.List[str]:
    # Read the mapping between the BUSCO lineages and their taxon_id
    with open(lineage_tax_ids) as file_in:
        lineage_tax_ids_dict: typing.Dict[int, str] = {}
        for line in file_in:
            arr = line.split()
            lineage_tax_ids_dict[int(arr[0])] = arr[1] + "_odb10"

    if requested_buscos:
        odb_arr = requested_buscos.split(",")
        valid_odbs = set(lineage_tax_ids_dict.values())
        for odb in odb_arr:
            if odb not in valid_odbs:
                print(f"Invalid BUSCO lineage: {odb}", file=sys.stderr)
                sys.exit(1)
    else:
        # Do the intersection to find the ancestors that have a BUSCO lineage
        odb_arr = [
            lineage_tax_ids_dict[anc_taxon_info.taxon_id]
            for anc_taxon_info in taxon_info.lineage
            if anc_taxon_info.taxon_id in lineage_tax_ids_dict
        ]

    return odb_arr


def get_assembly_info(accession: str) -> typing.Dict[str, typing.Union[str, int]]:
    response = requests.get(NCBI_DATASETS_API % accession).json()
    if response["total_count"] != 1:
        print(f"Assembly not found: {accession}", file=sys.stderr)
        sys.exit(1)
    assembly_report = response["reports"][0]
    assembly_info = assembly_report["assembly_info"]
    scaffold_count = assembly_report["assembly_stats"]["number_of_component_sequences"]
    scaffold_count += assembly_report["assembly_stats"].get("number_of_organelles", 0)
    span = int(assembly_report["assembly_stats"]["total_sequence_length"])
    span += sum(int(oi["total_seq_length"]) for oi in assembly_report.get("organelle_info", []))
    d = {
        "accession": accession,
        "alias": assembly_info["assembly_name"],
        "bioproject": assembly_info["bioproject_accession"],
        "level": assembly_info["assembly_level"].lower(),
        "scaffold-count": scaffold_count,
        "span": span,
    }
    if "biosample" in assembly_info:
        d["biosample"] = assembly_info["biosample"]["accession"]
    if "wgs_info" in assembly_report:
        d["prefix"] = assembly_report["wgs_info"]["wgs_project_accession"]
    return d


def get_sequence_report(accession: str):
    response = requests.get(NCBI_SEQUENCE_API % accession).json()
    if not response["reports"]:
        print(f"Assembly not found: {accession}", file=sys.stderr)
        sys.exit(1)
    sequence_report = response["reports"]
    for rec in sequence_report:
        if rec["role"] == "assembled-molecule":
            rec["assembly_level"] = rec["assigned_molecule_location_type"]
        else:
            rec["assembly_level"] = "na"
            rec["chr_name"] = "na"
    return sequence_report


def adjust_taxon_id(nt: str, taxon_info: TaxonInfo) -> int:
    con = sqlite3.connect(os.path.join(nt, "taxonomy4blast.sqlite3"))
    cur = con.cursor()
    for taxon_id in [taxon_info.taxon_id] + [anc_taxon_info.taxon_id for anc_taxon_info in taxon_info.lineage]:
        res = cur.execute("SELECT * FROM TaxidInfo WHERE taxid = ?", (taxon_id,))
        if res.fetchone():
            return taxon_id


def datatype_to_platform(s):
    if s == "ont":
        return "OXFORD_NANOPORE"
    elif s.startswith("pacbio"):
        return "PACBIO_SMRT"
    elif s in ["hic", "illumina"]:
        return "ILLUMINA"
    else:
        return "OTHER"


def print_yaml(
    file_out,
    assembly_info: typing.Dict[str, typing.Union[str, int]],
    taxon_info: TaxonInfo,
    classification: typing.Dict[str, str],
    reads,
    blastp,
    blastx,
    blastn,
    taxdump,
):
    data = {
        "assembly": assembly_info,
        "reads": {
            "paired": [],
            "single": [],
        },
        "revision": 1,
        "settings": {
            "blast_chunk": 100000,
            "blast_max_chunks": 10,
            "blast_min_length": 1000,
            "blast_overlap": 0,
            "stats_chunk": 1000,
            "stats_windows": [0.1, 0.01, 100000, 1000000],
            "taxdump": taxdump,
            "tmp": "/tmp",
        },
        "similarity": {
            "blastn": {
                "name": "nt",
                "path": blastn,
            },
            "defaults": {"evalue": 1e-10, "import_evalue": 1e-25, "max_target_seqs": 10, "taxrule": "buscogenes"},
            "diamond_blastp": {
                "import_max_target_seqs": 100000,
                "name": "reference_proteomes",
                "path": blastp,
                "taxrule": "blastp=buscogenes",
            },
            "diamond_blastx": {
                "name": "reference_proteomes",
                "path": blastx,
            },
        },
        "taxon": {
            "taxid": str(taxon_info.taxon_id),  # str on purpose, to mimic blobtoolkit
            "name": taxon_info.organism_name,
            **classification,
        },
        "version": 2,
    }

    for id, datatype, layout, path in reads:
        data["reads"][layout.lower()].append(
            {
                "prefix": id,
                "file": path,
                "plaform": datatype_to_platform(datatype),
            }
        )

    out_dir = os.path.dirname(file_out)
    make_dir(out_dir)

    with open(file_out, "w") as fout:
        yaml.dump(data, fout)


def print_csv(file_out, taxon_id: int, odb_arr: typing.List[str]):
    out_dir = os.path.dirname(file_out)
    make_dir(out_dir)

    with open(file_out, "w") as fout:
        print("taxon_id", taxon_id, sep=",", file=fout)
        for odb_val in odb_arr:
            print("busco_lineage", odb_val, sep=",", file=fout)


def print_tsvs(output_prefix, sequence_report):
    categories_tsv = f"{output_prefix}.categories.tsv"
    synonyms_tsv = f"{output_prefix}.synonyms.tsv"
    with open(categories_tsv, "w") as fhc:
        with open(synonyms_tsv, "w") as fhs:
            print("identifier", "assembly_role", "assembly_level", "assembly_unit", sep="\t", file=fhc)
            print("identifier", "name", "assigned_name", "refseq_accession", sep="\t", file=fhs)
            for rec in sequence_report:
                print(
                    rec["genbank_accession"],
                    rec["role"],
                    rec["assembly_level"],
                    rec["assembly_unit"],
                    sep="\t",
                    file=fhc,
                )
                print(
                    rec["genbank_accession"],
                    rec["sequence_name"],
                    rec["chr_name"],
                    rec.get("refseq_accession", "na"),
                    sep="\t",
                    file=fhs,
                )


def main(args=None):
    args = parse_args(args)

    assembly_info: typing.Dict[str, typing.Union[str, int]]
    if args.accession:
        assembly_info = get_assembly_info(args.accession)
        sequence_report = get_sequence_report(args.accession)
    else:
        assembly_info = {"level": "scaffold"}
        sequence_report = None
    assembly_info["file"] = args.fasta

    taxon_info = fetch_taxon_info(args.taxon_query)
    classification = get_classification(taxon_info)
    odb_arr = get_odb(taxon_info, args.lineage_tax_ids, args.busco)
    taxon_id = adjust_taxon_id(args.nt, taxon_info)

    if sequence_report:
        print_tsvs(args.output_prefix, sequence_report)

    reads = zip(args.read_id, args.read_type, args.read_layout, args.read_path)

    print_yaml(
        f"{args.output_prefix}.yaml",
        assembly_info,
        taxon_info,
        classification,
        reads,
        args.blastp,
        args.blastx,
        args.blastn,
        args.taxdump,
    )
    print_csv(f"{args.output_prefix}.csv", taxon_id, odb_arr)


if __name__ == "__main__":
    sys.exit(main())
