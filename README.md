# ![sanger-tol/blobtoolkit](docs/images/sanger-tol-blobtoolkit_logo.png)

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/sanger-tol/blobtoolkit)
[![GitHub Actions CI Status](https://github.com/sanger-tol/blobtoolkit/actions/workflows/nf-test.yml/badge.svg)](https://github.com/sanger-tol/blobtoolkit/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/sanger-tol/blobtoolkit/actions/workflows/linting.yml/badge.svg)](https://github.com/sanger-tol/blobtoolkit/actions/workflows/linting.yml)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.7949058-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.7949058)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.5.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.5.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/sanger-tol/blobtoolkit)

## Introduction

**sanger-tol/blobtoolkit** is a bioinformatics pipeline that can be used to identify and analyse non-target DNA for eukaryotic genomes.
It takes a samplesheet of BAM/CRAM/FASTQ/FASTA files as input, calculates genome statistics, coverage and completeness information, combines them in a TSV file by window size to create a BlobDir dataset and static plots.

1. Calculate genome statistics in windows ([`fastawindows`](https://github.com/tolkit/fasta_windows))
2. Calculate Coverage ([`blobtk/depth`](https://github.com/blobtoolkit/blobtk))
3. Determine the appropriate BUSCO lineages from the taxonomy.
4. Run BUSCO ([`busco`](https://busco.ezlab.org/))
5. Extract BUSCO genes ([`blobtoolkit/extractbuscos`](https://github.com/blobtoolkit/blobtoolkit))
6. Run Diamond BLASTp against extracted BUSCO genes ([`diamond/blastp`](https://github.com/bbuchfink/diamond))
7. Run BLASTx against sequences with no hit ([`diamond/blastx`](https://github.com/bbuchfink/diamond))
8. Run BLASTn against sequences still with no hit ([`blast/blastn`](https://www.ncbi.nlm.nih.gov/books/NBK131777/))
9. Count BUSCO genes ([`blobtoolkit/countbuscos`](https://github.com/blobtoolkit/blobtoolkit))
10. Generate combined sequence stats across various window sizes ([`blobtoolkit/windowstats`](https://github.com/blobtoolkit/blobtoolkit))
11. Import analysis results into a BlobDir dataset ([`blobtoolkit/blobdir`](https://github.com/blobtoolkit/blobtoolkit))
12. Create static plot images ([`blobtk/images`](https://github.com/blobtoolkit/blobtk))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,datatype,datafile,library_layout
mMelMel3_hic,hic,GCA_922984935.2.hic.mMelMel3.cram,PAIRED
mMelMel1,illumina,GCA_922984935.2.illumina.mMelMel1.cram,PAIRED
mMelMel3_ont,ont,GCA_922984935.2.ont.mMelMel3.cram,SINGLE
```

Each row represents a read set (aligned or not).
The first column (sample name) must be unique.
If you have multiple read sets from the same actual sample, make sure you edit the sample names to make them unique.
The datatype refers to the sequencing technology used to generate the underlying raw data and follows a controlled vocabulary (`ont`, `hic`, `pacbio`, `pacbio_clr`, `illumina`).
The library layout indicates whether the reads are paired or single.
The aligned read files can be generated using the [sanger-tol/readmapping](https://github.com/sanger-tol/readmapping) pipeline.

Now, you can run the pipeline using:

```bash
nextflow run sanger-tol/blobtoolkit \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --fasta genome.fasta \
   --accession GCA_XXXXXXXXX.X \
   --taxon XXXX \
   --taxdump /path/to/taxdump/database \
   --blastp /path/to/diamond/database \
   --blastn /path/to/blastn/nt.nal \
   --blastx /path/to/blastx/database
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details, please refer to the [usage documentation](https://pipelines.tol.sanger.ac.uk/blobtoolkit/usage) and the [parameter documentation](https://pipelines.tol.sanger.ac.uk/blobtoolkit/parameters).

## BLAST Database Configuration

### BLASTn Database Requirements

The `--blastn` parameter requires a direct path to a **BLAST database file**, either a `.nal` (alias) file or a `.nin` (index) file. The pipeline validates that all required companion files are present.

#### Supported File Types:

1. **`.nal` file (preferred)** - BLAST alias file:

   ```bash
   --blastn /path/to/databases/nt.nal
   ```

2. **`.nin` file (fallback)** - BLAST index file (when .nal is not available):

   ```bash
   --blastn /path/to/databases/nt.nin
   ```

3. **Compressed archive** (for CI/testing):
   ```bash
   --blastn https://example.com/path/to/nt_database.tar.gz
   ```

#### Database Structure Requirements:

##### When using a `.nal` file:

The directory must contain all companion files with the same prefix:

- `db_name.nal` (alias file - the file you point to)
- `db_name.nin` or `db_name.##.nin` (index file(s))
- `db_name.nhr` or `db_name.##.nhr` (header file(s))
- `db_name.nsq` or `db_name.##.nsq` (sequence file(s))

##### When using a `.nin` file:

The directory must contain companion files with the same prefix:

- `db_name.nin` or `db_name.##.nin` (index file - the file you point to)
- `db_name.nhr` or `db_name.##.nhr` (header file(s))
- `db_name.nsq` or `db_name.##.nsq` (sequence file(s))

**Note**: `##` represents numbers like `00`, `01`, `02`, etc. for large databases split into multiple files.

#### Example Directory Structures:

##### Single File Pattern:

```
/data/blast/nt/
├── nt.nal                   # Point --blastn here
├── nt.nin                   # Required companion files
├── nt.nhr
├── nt.nsq
├── taxdb.btd                # Optional taxonomy files
└── taxonomy4blast.sqlite3
```

##### Numbered File Pattern (Large Databases):

```
/data/blast/nt/
├── nt.nal                   # Point --blastn here
├── nt.00.nin                # Required numbered companion files
├── nt.00.nhr
├── nt.00.nsq
├── nt.01.nin
├── nt.01.nhr
├── nt.01.nsq
├── nt.02.nin
├── nt.02.nhr
├── nt.02.nsq
├── taxdb.btd
├── taxdb.bti
└── taxonomy4blast.sqlite3
```

##### Using .nin file (when .nal is not available):

```
/data/blast/nt/
├── nt.nin                   # Point --blastn here (no .nal file)
├── nt.nhr                   # Required companion files
├── nt.nsq
├── taxdb.btd
├── taxdb.bti
└── taxonomy4blast.sqlite3
```

#### Troubleshooting:

- **Error: "Invalid BLAST database path"** - Ensure you're pointing to either a `.nal` or `.nin` file, not a directory
- **Error: "Missing required files"** - Verify that all companion files (`.nin`, `.nhr`, `.nsq`) exist with the same prefix
- **Error: "BLAST database appears incomplete"** - Check that all required BLAST database components are present
- **Error: "File not found"** - Verify the file path is correct and the file exists

#### Migration from Previous Versions:

If you were previously using `--blastn /path/to/taxonomy4blast.sqlite3`, you now need to:

1. Use `--blastn /path/to/nt.nal` (if available), or
2. Use `--blastn /path/to/nt.nin` (if .nal is not available)

## Pipeline output

For more details about the output files and reports, please refer to the [output documentation](https://pipelines.tol.sanger.ac.uk/blobtoolkit/output).

## Credits

sanger-tol/blobtoolkit was written in Nextflow by:

- [Alexander Ramos Diaz](https://github.com/alxndrdiaz)
- [Zaynab Butt](https://github.com/zb32)
- [Priyanka Surana](https://github.com/priyanka-surana)
- [Matthieu Muffato](https://github.com/muffato)
- [Tyler Chafin](https://github.com/tkchafin)
- [Yumi Sims](https://github.com/yumisims)
- [Damon-Lee Bernard Pointon](https://github.com/DLBPointon)

The original design and coding for [BlobToolKit software and Snakemake pipeline](https://github.com/blobtoolkit/blobtoolkit) was done by [Richard Challis](https://github.com/rjchallis) and [Sujai Kumar](https://github.com/sujaikumar).

We thank the following people for their extensive assistance in the development of this pipeline:

- [Guoying Qi](https://github.com/gq1)
- [Bethan Yates](https://github.com/BethYates)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

If you use sanger-tol/blobtoolkit for your analysis, please cite it using the following DOI: [10.5281/zenodo.7949058](https://doi.org/10.5281/zenodo.7949058)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
