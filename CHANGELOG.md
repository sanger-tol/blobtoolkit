# sanger-tol/blobtoolkit: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [[1.0.0](https://github.com/sanger-tol/blobtoolkit/releases/tag/1.0.0)] – ??? - [2023-04-06]

### Enhancements & fixes

Initial release of sanger-tol/blobtoolkit :tada:

This release marks the point where the pipeline was moved from Snakemake at [blobtoolkit/blobtoolkit](https://github.com/blobtoolkit/blobtoolkit) over to Nextflow DSL2 at [sanger-tol/blobtoolkit](https://github.com/sanger-tol/blobtoolkit). There are two subworkflows in the Snakemake version that are still being implemented in Nextflow – `diamond_blastx` and `blastn`.

### Parameters

| Old parameter | New parameter   |
| ------------- | --------------- |
|               | --input         |
|               | --fasta         |
|               | --accession     |
|               | --taxon         |
|               | --taxa_file     |
|               | --yaml          |
|               | --blastp_outext |
|               | --blastp_cols   |
|               | --taxdump       |
|               | --busco         |
|               | --uniprot       |

> **NB:** Parameter has been **updated** if both old and new parameter information is present. </br> **NB:** Parameter has been **added** if just the new parameter information is present. </br> **NB:** Parameter has been **removed** if new parameter information isn't present.

### Software dependencies

Note, since the pipeline is using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference. Only `Docker` or `Singularity` containers are supported, `conda` is not supported.

| Dependency    | Old version | New version |
| ------------- | ----------- | ----------- |
| blobtoolkit   |             | 4.1.4       |
| busco         |             | 5.4.3       |
| fasta_windows |             | 0.2.4       |
| goat          |             | 0.2.0       |
| gunzip        |             | 1.10        |
| mosdepth      |             | 0.3.3       |
| nextflow      |             | 22.10.6     |
| python        |             | 3.10.6      |
| samtools      |             | 1.15.1      |
| tar           |             | 1.30        |
| yaml          |             | 6.0         |

> **NB:** Dependency has been **updated** if both old and new version information is present. </br> **NB:** Dependency has been **added** if just the new version information is present. </br> **NB:** Dependency has been **removed** if version information isn't present.
