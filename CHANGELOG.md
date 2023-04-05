# nf-core/blobtoolkit: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [[1.0.0](https://github.com/sanger-tol/blobtoolkit/releases/tag/1.0.0)] – ??? - [2023-04-06]

### Enhancements & fixes

Initial release of sanger-tol/blobtoolkit :tada:

+ Required input: (1) Tree of Life organism ID (ToLID, --input) or (2) FASTA genome file (--fasta) and `samplesheet.csv` (--input) containing paths to aligned reads in CRAM format, see description of *Full samplesheet* in `docs/usage.md`.
+ Lineages to run BUSCO are obtained using [goat_taxonsearch module](https://nf-co.re/modules/goat_taxonsearch) instead of reading them from `config.yaml` file (--taxon: NCBI taxonomy ID or species binomial name).
+ Versions of tools are collected using [custom_dumpsoftwareversions module](https://nf-co.re/modules/custom_dumpsoftwareversions) instead of reading them from `config.yaml` file.
+ The pipeline requires aligned reads, `minimap2` is no longer included.
+ A GCA accession ID or YAML file is required to collect metadata (--accession or --yaml).
+ [fastawindows module](https://nf-co.re/modules/fastawindows) computes statistics across windows over the genome FASTA file.
+ Parameters for [diamond_blastp module](https://nf-co.re/modules/diamond_blastp) should be specified in `nextflow.config` file (). 

This release marks the point where the pipeline was moved from Snakemake at [blobtoolkit/blobtoolkit](https://github.com/blobtoolkit/blobtoolkit) over to Nextflow DSL2 at [sanger-tol/blobtoolkit](https://github.com/sanger-tol/blobtoolkit).

### Parameters

| Old parameter | New parameter |
| ------------- | ------------- |
|               | --accession   |
|               | --blastp_outext |
|               | --blastp_cols   |
|               | --fasta       |
|               | --input       |
|               | --taxon       |
|               | --yaml        |

> **NB:** Parameter has been **updated** if both old and new parameter information is present. </br> **NB:** Parameter has been **added** if just the new parameter information
> is present. </br> **NB:** Parameter has been **removed** if new parameter information isn't present.

### Software dependencies

Note, since the pipeline is using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is e
ntirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been liste
d below for reference.

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
|blastn      | 2.12.0+     |             |
|blobtools   | 4.0.7       |4.1.2        |
|busco       | 5.3.2       |5.4.3        |
|diamond     | 2.0.15      |2.0.15       |
|fasta_windows |           | 0.2.4       |
|goat-cli    | 0.2.0       |0.2.0        |
|gunzip      |             |1.10         |
|minimap2    | 2.24-r1122  |             |
|mosdepth    | 0.3.3       |0.3.3        |
|multiqc     |             |1.13         |
|python      | 3.9.13      |3.9.13       |
|samtools    | 1.15.1      |1.15.1       |
|seqtk       | 1.3-r106    |             |

> **NB:** Dependency has been **updated** if both old and new version information is present. </br> **NB:** Dependency has been **added** if just the new version information is
> present. </br> **NB:** Dependency has been **removed** if version information isn't present.
