# ![sanger-tol/blobtoolkit](docs/images/sanger-tol-blobtoolkit_logo.png)

[![GitHub Actions CI Status](https://github.com/sanger-tol/blobtoolkit/workflows/nf-core%20CI/badge.svg)](https://github.com/sanger-tol/blobtoolkit/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/sanger-tol/blobtoolkit/workflows/nf-core%20linting/badge.svg)](https://github.com/sanger-tol/blobtoolkit/actions?query=workflow%3A%22nf-core+linting%22)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.7949058-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.7949058)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/sanger-tol/blobtoolkit)

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
8. Run BLASTn against sequences still with not hit ([`blast/blastn`](https://www.ncbi.nlm.nih.gov/books/NBK131777/))
9. Count BUSCO genes ([`blobtoolkit/countbuscos`](https://github.com/blobtoolkit/blobtoolkit))
10. Generate combined sequence stats across various window sizes ([`blobtoolkit/windowstats`](https://github.com/blobtoolkit/blobtoolkit))
11. Imports analysis results into a BlobDir dataset ([`blobtoolkit/blobdir`](https://github.com/blobtoolkit/blobtoolkit))
12. Create static plot images ([`blobtk/images`](https://github.com/blobtoolkit/blobtk))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,datatype,datafile,library_layout
mMelMel3,hic,GCA_922984935.2.hic.mMelMel3.cram,PAIRED
mMelMel1,illumina,GCA_922984935.2.illumina.mMelMel1.cram,PAIRED
mMelMel3,ont,GCA_922984935.2.ont.mMelMel3.cram,SINGLE
```

Each row represents an aligned file.
Rows with the same sample identifier are considered technical replicates.
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
   --blastn /path/to/blastn/database \
   --blastx /path/to/blastx/database
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

For more details, please refer to the [usage documentation](https://pipelines.tol.sanger.ac.uk/blobtoolkit/usage) and the [parameter documentation](https://pipelines.tol.sanger.ac.uk/blobtoolkit/parameters).

## Pipeline output

For more details about the output files and reports, please refer to the [output documentation](https://pipelines.tol.sanger.ac.uk/blobtoolkit/output).

## Credits

sanger-tol/blobtoolkit was written in Nextflow by [Alexander Ramos Diaz](https://github.com/alxndrdiaz), [Zaynab Butt](https://github.com/zb32), [Matthieu Muffato](https://github.com/muffato), and [Priyanka Surana](https://github.com/priyanka-surana). The orignal design and coding for [BlobToolKit software and Snakemake pipeline](https://github.com/blobtoolkit/blobtoolkit) was done by [Richard Challis](https://github.com/rjchallis) and [Sujai Kumar](https://github.com/sujaikumar).

We thank the following people for their assistance in the development of this pipeline:

- [Guoying Qi](https://github.com/gq1)
- [Bethan Yates](https://github.com/BethYates)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

If you use sanger-tol/blobtoolkit for your analysis, please cite it using the following doi: [10.5281/zenodo.7949058](https://doi.org/10.5281/zenodo.7949058)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
