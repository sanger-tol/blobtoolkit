# sanger-tol/blobtoolkit: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [[1.0.0](https://github.com/sanger-tol/blobtoolkit/releases/tag/1.0.0)] – Psyduck – [2024-XX-YY]

The pipeline is now considered to be a complete and suitable replacement for the Snakemake version.

- Fetch information about the chromosomes of the assemblies. Used to power
  "grid plots".
- Fill in accurate read information in the blobDir. Users are now reqiured
  to indicate whether the reads are paired or single.

## [[0.6.0](https://github.com/sanger-tol/blobtoolkit/releases/tag/0.6.0)] – Bellsprout – [2024-09-13]

The pipeline has now been validated for draft (unpublished) assemblies.

- The pipeline now queries the NCBI database instead of GoaT to establish the
  taxonomic classification of the species and the relevant Busco lineages.
  In case the taxon_id is not found, the pipeline falls back to GoaT, which
  is aware of upcoming taxon_ids in ENA.
- New `--busco_lineages` parameter to choose specific Busco lineages instead of
  automatically selecting based on the taxonomy.
- All parameters are now passed the regular Nextflow way. There is no support
  for the original Yaml configuration files of the Snakemake version.
- New option `--skip_taxon_filtering` to skip the taxon filtering in blast searches.
  Mostly relevant for draft assemblies.
- Introduced the `--use_work_dir_as_temp` parameter to avoid leaving files in `/tmp`.

### Parameters

| Old parameter | New parameter          |
| ------------- | ---------------------- |
| --yaml        |                        |
|               | --busco_lineages       |
|               | --skip_taxon_filtering |
|               | --use_work_dir_as_temp |

> **NB:** Parameter has been **updated** if both old and new parameter information is present. </br> **NB:** Parameter has been **added** if just the new parameter information is present. </br> **NB:** Parameter has been **removed** if new parameter information isn't present.

### Software dependencies

Note, since the pipeline is using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference. Only `Docker` or `Singularity` containers are supported, `conda` is not supported.

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| goat       | 0.2.5       |             |

## [[0.5.1](https://github.com/sanger-tol/blobtoolkit/releases/tag/0.5.1)] – Snorlax (patch 1) – [2024-08-22]

### Enhancements & fixes

- Bugfix: skip BLASTN if there are no chunks to align

## [[0.5.0](https://github.com/sanger-tol/blobtoolkit/releases/tag/0.5.0)] – Snorlax – [2024-07-31]

General tidy up of the configuration and the pipeline

### Enhancements & fixes

- Increased the resources for blastn
- Removed some options that were not used or not needed
- All relevant outputs are now copied to the output directory
- Fixed some blast parameters to match the behaviour of the Snakemake pipeline
- Fixed parsing of samplesheets from fetchngs to capture correct data type

### Parameters

| Old parameter   | New parameter |
| --------------- | ------------- |
| --taxa_file     |               |
| --blastp_outext |               |
| --blastp_cols   |               |
| --blastx_outext |               |
| --blastx_cols   |               |

> **NB:** Parameter has been **updated** if both old and new parameter information is present. </br> **NB:** Parameter has been **added** if just the new parameter information is present. </br> **NB:** Parameter has been **removed** if new parameter information isn't present.

## [[0.4.0](https://github.com/sanger-tol/blobtoolkit/releases/tag/0.4.0)] – Buneary – [2024-04-17]

The pipeline has now been validated on dozens of genomes, up to 11 Gbp.

### Enhancements & fixes

- Upgraded the version of `blobtools`, which enables a better reporting of
  wrong accession numbers and a better handling of oddities in input files.
- Files in the output blobdir are now compressed.
- All modules handling blobdirs can now be cached.
- Large genomes supported, up to at least 11 Gbp.
- Allow all variations of FASTA and FASTQ extensions for input.
- More fields included in the trace files.
- All nf-core modules updated

### Software dependencies

Note, since the pipeline is using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference. Only `Docker` or `Singularity` containers are supported, `conda` is not supported.

| Dependency  | Old version   | New version       |
| ----------- | ------------- | ----------------- |
| blobtoolkit | 4.3.3         | 4.3.9             |
| blast       | 2.14.0        | 2.15.0 and 2.14.1 |
| multiqc     | 1.17 and 1.18 | 1.20 and 1.21     |
| samtools    | 1.18          | 1.19.2            |
| seqtk       | 1.3           | 1.4               |

> **NB:** Dependency has been **updated** if both old and new version information is present. </br> **NB:** Dependency has been **added** if just the new version information is present. </br> **NB:** Dependency has been **removed** if version information isn't present.

## [[0.3.0](https://github.com/sanger-tol/blobtoolkit/releases/tag/0.3.0)] – Poliwag – [2024-02-09]

The pipeline has now been validated on five genomes, all under 100 Mbp: a
sponge, a platyhelminth, and three fungi.

### Enhancements & fixes

- Fixed the conditional runs of blastn
- Fixed the generation of the no-hit list
- Fixed the conversion of the unaligned input files to Fasta
- Fixed the documentation about preparing the NT database
- Fixed the detection of the NT database in the nf-core module
- The pipeline now supports samplesheets generated by the
  [nf-core/fetchngs](https://nf-co.re/fetchngs) pipeline by passing the
  `--fetchngs_samplesheet true` option.
- FastQ files can bypass the conversion to Fasta
- Fixed missing BUSCO results from the blobdir (only 1 BUSCO was loaded)
- Fixed the default category used to colour the blob plots
- Fixed the output directory of the images
- Added an option to select the format of the images (PNG or SVG)

### Parameters

| Old parameter | New parameter          |
| ------------- | ---------------------- |
|               | --fetchngs_samplesheet |
|               | --image_format         |

> **NB:** Parameter has been **updated** if both old and new parameter information is present. </br> **NB:** Parameter has been **added** if just the new parameter information is present. </br> **NB:** Parameter has been **removed** if new parameter information isn't present.

### Software dependencies

Note, since the pipeline is using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference. Only `Docker` or `Singularity` containers are supported, `conda` is not supported.

| Dependency  | Old version | New version |
| ----------- | ----------- | ----------- |
| blobtoolkit | 4.3.2       | 4.3.3       |

> **NB:** Dependency has been **updated** if both old and new version information is present. </br> **NB:** Dependency has been **added** if just the new version information is present. </br> **NB:** Dependency has been **removed** if version information isn't present.

## [[0.2.0](https://github.com/sanger-tol/blobtoolkit/releases/tag/0.2.0)] – Pikachu – [2023-12-22]

### Enhancements & fixes

- Template updated to nf-core/tools 2.11.1
- Includes all subworkflows in the [Snakemake version](https://github.com/blobtoolkit/blobtoolkit)
- Added blastx and blastn subworkflows
- Replaced mosdepth with blobtk depth
- Updated config creation script

### Parameters

| Old parameter | New parameter   |
| ------------- | --------------- |
|               | --mask          |
|               | --align         |
| --uniprot     | --blastp        |
|               | --blastx        |
|               | --blastn        |
|               | --blastx_outext |
|               | --blastx_cols   |

> **NB:** Parameter has been **updated** if both old and new parameter information is present. </br> **NB:** Parameter has been **added** if just the new parameter information is present. </br> **NB:** Parameter has been **removed** if new parameter information isn't present.

### Software dependencies

Note, since the pipeline is using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference. Only `Docker` or `Singularity` containers are supported, `conda` is not supported.

| Dependency   | Old version | New version |
| ------------ | ----------- | ----------- |
| blobtoolkit  | 4.1.4       | 4.3.2       |
| busco        | 5.4.3       | 5.5.0       |
| goat         | 0.2.0       | 0.2.5       |
| mosdepth     | 0.3.3       |             |
| nextflow     | 22.10.6     | 23.10.0     |
| python       | 3.10.6      | 3.12.0      |
| samtools     | 1.15.1      | 1.18        |
| tar          | 1.30        |             |
| yaml         | 6.0         | 6.0.1       |
| blobtk       | 0.3.3       | 0.5.1       |
| diamond      | 2.0.15      | 2.1.8       |
| minimap2     |             | 2.24-r1122  |
| blast        |             | 2.14.1      |
| windowmasker |             | 2.14.0      |

> **NB:** Dependency has been **updated** if both old and new version information is present. </br> **NB:** Dependency has been **added** if just the new version information is present. </br> **NB:** Dependency has been **removed** if version information isn't present.

## [[0.1.0](https://github.com/sanger-tol/blobtoolkit/releases/tag/0.1.0)] – Vaporeon – [2023-05-18]

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
