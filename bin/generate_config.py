#!/usr/bin/env python3

import argparse
import dataclasses
import os
import sys
import requests
import typing
import yaml

NCBI_TAXONOMY_API = "https://api.ncbi.nlm.nih.gov/datasets/v1/taxonomy/taxon/%s"

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
    parser.add_argument("FASTA", help="Path to the Fasta file of the assembly.")
    parser.add_argument("TAXON_QUERY", help="Query string/integer for this taxon.")
    parser.add_argument("LINEAGE_TAX_IDS", help="Mapping between BUSCO lineages and taxon IDs.")
    parser.add_argument("YML_OUT", help="Output YML file.")
    parser.add_argument("CSV_OUT", help="Output CSV file.")
    parser.add_argument("--busco", dest="REQUESTED_BUSCOS", help="Requested BUSCO lineages.", default=None)
    parser.add_argument("--version", action="version", version="%(prog)s 1.1")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        os.makedirs(path, exist_ok=True)


@dataclasses.dataclass
class TaxonInfo:
    taxon_id: int
    organism_name: str
    rank: typing.Optional[str]
    lineage: typing.List[int]


def make_taxon_info(taxon_name: typing.Union[str, int]) -> TaxonInfo:
    # Using API, get the taxon_ids of the species and all parents
    response = requests.get(NCBI_TAXONOMY_API % taxon_name).json()
    body = response["taxonomy_nodes"][0]["taxonomy"]
    return TaxonInfo(body["tax_id"], body["organism_name"], body.get("rank"), body["lineage"])


def get_classification(taxon_info: TaxonInfo) -> typing.Dict[str, str]:
    ancestors: typing.Dict[str, str] = {}
    for anc_taxon_id in taxon_info.lineage[1:]:  # skip the root taxon_id=1
        anc_taxon_info = make_taxon_info(anc_taxon_id)
        if anc_taxon_info.rank:
            ancestors[anc_taxon_info.rank.lower()] = anc_taxon_info.organism_name
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
            lineage_tax_ids_dict[taxon_id]
            for taxon_id in reversed(taxon_info.lineage)
            if taxon_id in lineage_tax_ids_dict
        ]

    return odb_arr


def print_yaml(
    file_out, fasta: str, taxon_info: TaxonInfo, classification: typing.Dict[str, str], odb_arr: typing.List[str]
):

    data = {
        "assembly": {
            # Other attributes can't be filled in for draft assemblies
            "file": fasta,
            "level": "scaffold",
        },
        "busco": {
            "basal_lineages": BUSCO_BASAL_LINEAGES,
            # "download_dir": <completely skipped because missing from final meta.json>
            "lineages": odb_arr + [lin for lin in BUSCO_BASAL_LINEAGES if lin not in odb_arr],
        },
        # TODO: reads
        "revision": 1,
        "settings": {
            # Only settings.stats_windows is mandatory, everything else is superfluous
            "blast_chunk": 100000,
            "blast_max_chunks": 10,
            "blast_min_length": 1000,
            "blast_overlap": 0,
            "stats_chunk": 1000,
            "stats_windows": [0.1, 0.01, 100000, 1000000],
            # "taxdump": <added by update_versions.py at the end>,
            "tmp": "/tmp",
        },
        "similarity": {
            # Only the presence similarity.diamond_blastx seems mandatory, everything else is superfluous
            "blastn": {
                "name": "nt",
                # "path": <added by update_versions.py at the end>,
            },
            "defaults": {"evalue": 1e-10, "import_evalue": 1e-25, "max_target_seqs": 10, "taxrule": "buscogenes"},
            "diamond_blastp": {
                "import_max_target_seqs": 100000,
                "name": "reference_proteomes",
                # "path": <added by update_versions.py at the end>,
                "taxrule": "blastp=buscogenes",
            },
            "diamond_blastx": {
                "name": "reference_proteomes",
                # "path": <added by update_versions.py at the end>,
            },
        },
        "taxon": {
            "taxid": str(taxon_info.taxon_id),  # str on purpose, to mimic blobtoolkit
            "name": taxon_info.organism_name,
            **classification,
        },
        "version": 2,
    }
    out_dir = os.path.dirname(file_out)
    make_dir(out_dir)

    with open(file_out, "w") as fout:
        yaml.dump(data, fout)


def print_csv(file_out, taxon_info: TaxonInfo, odb_arr: typing.List[str]):

    out_dir = os.path.dirname(file_out)
    make_dir(out_dir)

    with open(file_out, "w") as fout:
        print("taxon_id", taxon_info.taxon_id, sep=",", file=fout)
        for odb_val in odb_arr:
            print("busco_lineage", odb_val, sep=",", file=fout)


def main(args=None):
    args = parse_args(args)

    taxon_info = make_taxon_info(args.TAXON_QUERY)
    classification = get_classification(taxon_info)
    odb_arr = get_odb(taxon_info, args.LINEAGE_TAX_IDS, args.REQUESTED_BUSCOS)

    print_yaml(args.YML_OUT, args.FASTA, taxon_info, classification, odb_arr)
    print_csv(args.CSV_OUT, taxon_info, odb_arr)


if __name__ == "__main__":
    sys.exit(main())
