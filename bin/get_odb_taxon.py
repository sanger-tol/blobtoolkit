#!/usr/bin/env python3

import argparse
import os
import sys
import requests


NCBI_TAXONOMY_API = "https://api.ncbi.nlm.nih.gov/datasets/v1/taxonomy/taxon/%s"


def parse_args(args=None):
    Description = "Get Taxon ID and ODB database value using NCBI API and BUSCO configuration file"

    parser = argparse.ArgumentParser(description=Description)
    parser.add_argument("TAXON_QUERY", help="Query string/integer for this taxon.")
    parser.add_argument("LINEAGE_TAX_IDS", help="Mapping between BUSCO lineages and taxon IDs.")
    parser.add_argument("FILE_OUT", help="Output CSV file.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        os.makedirs(path, exist_ok=True)


def get_odb(taxon_query, lineage_tax_ids, file_out):
    # Read the mapping between the BUSCO lineages and their taxon_id
    with open(lineage_tax_ids) as file_in:
        lineage_tax_ids_dict = {}
        for line in file_in:
            arr = line.split()
            lineage_tax_ids_dict[int(arr[0])] = arr[1]

    # Using API, get the taxon_ids of the species and all parents
    response = requests.get(NCBI_TAXONOMY_API % taxon_query).json()
    this_taxon_id = response["taxonomy_nodes"][0]["taxonomy"]["tax_id"]
    ancestor_taxon_ids = response["taxonomy_nodes"][0]["taxonomy"]["lineage"]

    # Do the intersection to find the ancestors that have a BUSCO lineage
    odb_arr = [
        lineage_tax_ids_dict[taxon_id] for taxon_id in reversed(ancestor_taxon_ids) if taxon_id in lineage_tax_ids_dict
    ]

    out_dir = os.path.dirname(file_out)
    make_dir(out_dir)

    with open(file_out, "w") as fout:
        print("taxon_id", this_taxon_id, sep=",", file=fout)
        for odb_val in odb_arr:
            print("busco_lineage", odb_val + "_odb10", sep=",", file=fout)


def main(args=None):
    args = parse_args(args)
    get_odb(args.TAXON_QUERY, args.LINEAGE_TAX_IDS, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
