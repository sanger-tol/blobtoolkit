#!/usr/bin/env python3

"""Script to update the software versions in meta.json"""

import argparse
import sys
import json
import yaml


def parse_args(args=None):
    Description = "Combine BED files to create window stats input file."

    parser = argparse.ArgumentParser(description=Description)
    parser.add_argument("--meta_in", help="Input JSON file.", required=True)
    parser.add_argument("--meta_out", help="Output JSON file.", required=True)
    parser.add_argument("--software", help="Input YAML file.", required=True)
    parser.add_argument("--blastp", help="Path to the blastp database", required=True)
    parser.add_argument("--blastx", help="Path to the blastx database", required=True)
    parser.add_argument("--blastn", help="Path to the blastn database", required=True)
    parser.add_argument("--taxdump", help="Path to the taxonomy database", required=True)
    parser.add_argument("--version", action="version", version="%(prog)s 1.2.0")
    return parser.parse_args(args)


def datatype_to_platform(s):
    if s == "ont":
        return "OXFORD_NANOPORE"
    elif s.startswith("pacbio"):
        return "PACBIO_SMRT"
    elif s in ["hic", "illumina"]:
        return "ILLUMINA"
    else:
        return "OTHER"


def update_meta(args):
    with open(args.meta_in) as fh:
        infile = json.load(fh)

    with open(args.software) as fh:
        versions = yaml.safe_load(fh)

    new_dict = dict()
    for k, v in versions.items():
        new_dict.update(v)

    infile["settings"]["pipeline"] = "https://github.com/sanger-tol/blobtoolkit"
    infile["settings"]["release"] = new_dict["sanger-tol/blobtoolkit"]

    del new_dict["sanger-tol/blobtoolkit"]
    infile["settings"]["software_versions"] = new_dict

    infile["settings"]["taxdump"] = args.taxdump
    for k in ["blastn", "diamond_blastp", "diamond_blastx"]:
        infile["similarity"].setdefault(k, {})
    infile["similarity"]["blastn"]["path"] = args.blastn
    infile["similarity"]["diamond_blastp"]["path"] = args.blastp
    infile["similarity"]["diamond_blastx"]["path"] = args.blastx
    return infile


def main(args=None):
    args = parse_args(args)

    data = update_meta(args)
    with open(args.meta_out, "w") as fh:
        json.dump(data, fh)


if __name__ == "__main__":
    sys.exit(main())
