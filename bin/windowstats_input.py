#!/usr/bin/env python3

import argparse
import os
import sys
import pandas as pd


def parse_args(args=None):
    Description = "Combine BED files to create window stats input file."

    parser = argparse.ArgumentParser(description=Description)
    parser.add_argument(
        "--freq", help="Frequence fasta windows input file", required=True
    )
    parser.add_argument(
        "--mononuc", help="Mononucleotide fasta windows input file", required=True
    )
    parser.add_argument(
        "--mosdepth", help="Mosdepth coverage input file", nargs="+", required=True
    )
    parser.add_argument(
        "--countbusco", help="BUSCO gene counts by region", required=True
    )
    parser.add_argument("--output", help="Output TSV file.", required=True)
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        os.makedirs(path, exist_ok=True)


def merge_all(freq, mononuc, mosdepth, countbusco):
    freq_fw = pd.read_csv(freq, sep="\t")
    mononuc_fw = pd.read_csv(mononuc, sep="\t")
    combo_fw = freq_fw.merge(mononuc_fw).rename(
        columns={"ID": "sequence", "GC_prop": "gc", "Prop_Ns": "n", "N": "ncount"}
    )

    count_df = pd.read_csv(countbusco, sep="\t").rename(columns={"ID": "sequence"})
    for f in mosdepth:
        tag = os.path.basename(f).replace(".regions.bed.gz", "")
        cov_df = pd.read_csv(
            f,
            compression="gzip",
            sep="\t",
            names=["sequence", "start", "end", tag + "_cov"],
        )
        count_df = count_df.merge(cov_df)

    combo_all = combo_fw.merge(count_df)
    return combo_all


def main(args=None):
    args = parse_args(args)

    out_dir = os.path.dirname(args.output)
    make_dir(out_dir)

    merge_all(args.freq, args.mononuc, args.mosdepth, args.countbusco).to_csv(
        args.output, sep="\t", index=False
    )


if __name__ == "__main__":
    sys.exit(main())
