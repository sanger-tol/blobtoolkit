#!/usr/bin/env python3

import argparse
import sys

from blobtools.lib.file_io import write_file
from blobtools.lib.taxdump import Taxdump


def parse_args(args=None):
    Description = "Parse and digest the taxdump files into a JSON structure, printed on stdout."

    parser = argparse.ArgumentParser(description=Description)
    parser.add_argument("taxdump", help="Path to the taxonomy database")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0")
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)

    taxdump = Taxdump(args.taxdump)
    write_file("STDOUT", taxdump.values_to_dict())


if __name__ == "__main__":
    sys.exit(main())
