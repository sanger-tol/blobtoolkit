#!/usr/bin/env python


"""Provide a command line tool to validate and transform tabular samplesheets."""


import argparse
import csv
import logging
import sys
from collections import Counter
from pathlib import Path

logger = logging.getLogger()


class RowChecker:
    """
    Define a service that can validate and transform each given row.

    Attributes:
        modified (list): A list of dicts, where each dict corresponds to a previously
            validated and transformed row. The order of rows is maintained.

    """

    VALID_FORMATS = (
        ".fastq.gz",
    )

    def __init__(
        self,
        accession_col="run_accession",
        model_col="instrument_model",
        platform_col="instrument_platform",
        library_col="library_strategy",
        file1_col="fastq_1",
        file2_col="fastq_2",
        **kwargs,
    ):
        """
        Initialize the row checker with the expected column names.

        Args:
            accession_col (str): The name of the column that contains the accession name
                (default "run_accession").
            model_col (str): The name of the column that contains the model name
                of the instrument (default "instrument_model").
            platform_col (str): The name of the column that contains the platform name
                of the instrument (default "instrument_platform").
            library_col (str): The name of the column that contains the strategy of the
                preparation of the library (default "library_strategy").
            file2_col (str): The name of the column that contains the second file path
                for the paired-end read data (default "fastq_2").
        """
        super().__init__(**kwargs)
        self._accession_col = accession_col
        self._model_col = model_col
        self._platform_col = platform_col
        self._library_col = library_col
        self._file1_col = file1_col
        self._file2_col = file2_col
        self._seen = set()
        self.modified = []

    def validate_and_transform(self, row):
        """
        Perform all validations on the given row.

        Args:
            row (dict): A mapping from column headers (keys) to elements of that row
                (values).

        """
        self._validate_accession(row)
        self._validate_file(row)
        self._seen.add((row[self._accession_col], row[self._file1_col]))
        self.modified.append(row)

    def _validate_accession(self, row):
        """Assert that the run accession name exists."""
        if len(row[self._accession_col]) <= 0:
            raise AssertionError("Run accession is required.")

    def _validate_file(self, row):
        """Assert that the datafile is non-empty and has the right format."""
        if len(row[self._file1_col]) <= 0:
            raise AssertionError("Data file is required.")
        self._validate_data_format(row[self._file1_col])
        if row[self._file2_col]:
            self._validate_data_format(row[self._file2_col])

    def _validate_data_format(self, filename):
        """Assert that a given filename has one of the expected FASTQ extensions."""
        if not any(filename.endswith(extension) for extension in self.VALID_FORMATS):
            raise AssertionError(
                f"The data file has an unrecognized extension: {filename}\n"
                f"It should be one of: {', '.join(self.VALID_FORMATS)}"
            )

    def validate_unique_accessions(self):
        """
        Assert that the combination of accession name and aligned filename is unique.

        In addition to the validation, also rename all accessions to have a suffix of _T{n}, where n is the
        number of times the same accession exist, but with different FASTQ files, e.g., multiple runs per experiment.

        """
        if len(self._seen) != len(self.modified):
            raise AssertionError("The pair of accession and file name must be unique.")
        seen = Counter()
        for row in self.modified:
            accession = row[self._accession_col]
            seen[accession] += 1
            row[self._accession_col] = f"{accession}_T{seen[accession]}"


def read_head(handle, num_lines=10):
    """Read the specified number of lines from the current position in the file."""
    lines = []
    for idx, line in enumerate(handle):
        if idx == num_lines:
            break
        lines.append(line)
    return "".join(lines)


def sniff_format(handle):
    """
    Detect the tabular format.

    Args:
        handle (text file): A handle to a `text file`_ object. The read position is
        expected to be at the beginning (index 0).

    Returns:
        csv.Dialect: The detected tabular format.

    .. _text file:
        https://docs.python.org/3/glossary.html#term-text-file

    """
    peek = read_head(handle)
    handle.seek(0)
    sniffer = csv.Sniffer()
    dialect = sniffer.sniff(peek)
    return dialect


def check_samplesheet(file_in, file_out):
    """
    Check that the tabular samplesheet has the structure expected by sanger-tol pipelines.

    Validate the general shape of the table, expected columns, and each row. Also add
    Args:
        file_in (pathlib.Path): The given tabular samplesheet. The format can be either
            CSV, TSV, or any other format automatically recognized by ``csv.Sniffer``.
        file_out (pathlib.Path): Where the validated and transformed samplesheet should
            be created; always in CSV format.

    Example:
        This function checks that the samplesheet follows the following structure,
        see also the `blobtoolkit samplesheet`_::

        sample,datatype,datafile
        sample1,hic,/path/to/file1.cram
        sample1,pacbio,/path/to/file2.cram
        sample1,ont,/path/to/file3.cram

    .. _blobtoolkit samplesheet:
        https://raw.githubusercontent.com/sanger-tol/blobtoolkit/main/assets/test/samplesheet.csv

    """
    required_columns = {"run_accession", "instrument_model", "instrument_platform", "library_strategy", "fastq_1", "fastq_2"}
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_in.open(newline="") as in_handle:
        reader = csv.DictReader(in_handle, dialect=sniff_format(in_handle))
        # Validate the existence of the expected header columns.
        if not required_columns.issubset(reader.fieldnames):
            req_cols = ", ".join(required_columns)
            logger.critical(f"The sample sheet **must** contain these column headers: {req_cols}.")
            sys.exit(1)
        # Validate each row.
        checker = RowChecker()
        for i, row in enumerate(reader):
            try:
                checker.validate_and_transform(row)
            except AssertionError as error:
                logger.critical(f"{str(error)} On line {i + 2}.")
                sys.exit(1)
        checker.validate_unique_accessions()
    header = list(reader.fieldnames)
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_out.open(mode="w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, header, delimiter=",")
        writer.writeheader()
        for row in checker.modified:
            writer.writerow(row)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Tabular input samplesheet in CSV or TSV format.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Transformed output samplesheet in CSV format.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s 1.0.0",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    check_samplesheet(args.file_in, args.file_out)


if __name__ == "__main__":
    sys.exit(main())
