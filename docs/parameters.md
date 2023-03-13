# Parameters

BlobToolKit Nextflow Pipeline.

## Input/output options

Define where the pipeline should find input data and save output data.

| Parameter       | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                | Type     | Default | Required | Hidden |
| --------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------- | ------- | -------- | ------ |
| `input`         | Path to comma-separated file containing information about the samples in the experiment. <details><summary>Help</summary><small>Provide a Tree of Life organism ID. Otherwise, you will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row. See (https://nf-co.re/variantcalling/usage#samplesheet-input).</small></details> | `string` |         |          |        |
| `outdir`        | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.                                                                                                                                                                                                                                                                                                                                                                                   | `string` |         |          |        |
| `email`         | Email address for completion summary. <details><summary>Help</summary><small>Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.</small></details>                                                                                                                                                | `string` |         |          |        |
| `multiqc_title` | MultiQC report title. Printed as page header, used for filename if not otherwise specified.                                                                                                                                                                                                                                                                                                                                                                                                                | `string` |         |          |        |

Additionaly a path to a YAML file or an accesion code (usually a GCA or draft identifier) should be provided through `params.yaml` or `params.accesion` (only one of them should be specified). Here is an example of how a YAML file should look like, all information can be obtained from [NCBI: Browse taxonomy](https://www.ncbi.nlm.nih.gov/data-hub/taxonomy/9662/) :

```
assembly:
  accession: GCA_922984935.2
  level: chromosome
  prefix: CAKLPM02
  scaffold-count: 538
  span: 2738694574
revision: 1
settings:
  software_versions:
    blastn: 2.12.0+
    blobtools: 4.0.7
    busco: 5.3.2
    diamond: 2.0.15
    minimap2: 2.24-r1122
    mosdepth: 0.3.3
    python: 3.9.13
    samtools: 1.15.1
    seqtk: 1.3-r106
    snakemake: 7.19.1
  stats_windows:
    - 0.1
    - 0.01
    - 100000
    - 1000000
similarity:
  diamond_blastp:
    evalue: 1.0e-10
    import_evalue: 1.0e-25
    import_max_target_seqs: 100000
    max_target_seqs: 10
    path: /blobtoolkit/databases/uniprot_2021_06
    taxrule: blastp=buscogenes
    name: reference_proteomes
  diamond_blastx:
    evalue: 1.0e-10
    import_evalue: 1.0e-25
    max_target_seqs: 10
    name: reference_proteomes
    path: /blobtoolkit/databases/uniprot_2021_06
    taxrule: buscogenes
taxon:
  class: Mammalia
  family: Mustelidae
  genus: Meles
  kingdom: Metazoa
  name: Meles meles
  order: Carnivora
  phylum: Chordata
  superkingdom: Eukaryota
  taxid: '9662'
version: 2
```

Parameters in stats_windows, diamond_blastp, and diamond_blastx are ignored and are kept in this YAML file only to allow compatibility with the blobltools subworkflow in the previous [blobtoolkit-pipeline](https://github.com/blobtoolkit/blobtoolkit/tree/main/src/blobtoolkit-pipeline/src) implementation. If you need to specify new values for these parameters, before running the pipeline, you can simply edit the `conf/modules.config` file.

## Reference genome options

Reference genome related files and options required for the workflow.

| Parameter         | Description                                                                                                                                                                                                                                                                                                                                                                                     | Type      | Default                    | Required | Hidden |
| ----------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------- | -------------------------- | -------- | ------ |
| `genome`          | Name of iGenomes reference. <details><summary>Help</summary><small>If using a reference genome configured in the pipeline using iGenomes, use this parameter to give the ID for the reference. This is then used to build the full paths for all required reference genome files e.g. `--genome GRCh38`. See the (https://nf-co.re/usage/reference_genomes) for more details.</small></details> | `string`  |                            |          |        |
| `fasta`           | Path to FASTA genome file. <details><summary>Help</summary><small>This parameter is _mandatory_ if `--genome` is not specified. If you don't have a BWA index available this will be generated for you automatically. Combine with `--save_reference` to save BWA index for future runs.</small></details>                                                                                      | `string`  |                            |          |        |
| `igenomes_base`   | Directory / URL base for iGenomes references.                                                                                                                                                                                                                                                                                                                                                   | `string`  | s3://ngi-igenomes/igenomes |          | True   |
| `igenomes_ignore` | Do not load the iGenomes reference config. <details><summary>Help</summary><small>Do not load `igenomes.config` when running the pipeline. You may choose this option if you observe clashes between custom parameters and those supplied in `igenomes.config`.</small></details>                                                                                                               | `boolean` |                            |          | True   |

## Institutional config options

Parameters used to describe centralised config profiles. These should not be edited.

| Parameter                    | Description                                                                                                                                                                                                                                                                                                                                                                                       | Type     | Default                                                  | Required | Hidden |
| ---------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------- | -------------------------------------------------------- | -------- | ------ |
| `custom_config_version`      | Git commit id for Institutional configs.                                                                                                                                                                                                                                                                                                                                                          | `string` | master                                                   |          | True   |
| `custom_config_base`         | Base directory for Institutional configs. <details><summary>Help</summary><small>If you're running offline, Nextflow will not be able to fetch the institutional config files from the internet. If you don't need them, then this is not a problem. If you do need them, you should download the files from the repo and tell Nextflow where to find them with this parameter.</small></details> | `string` | https://raw.githubusercontent.com/nf-core/configs/master |          | True   |
| `config_profile_name`        | Institutional config name.                                                                                                                                                                                                                                                                                                                                                                        | `string` |                                                          |          | True   |
| `config_profile_description` | Institutional config description.                                                                                                                                                                                                                                                                                                                                                                 | `string` |                                                          |          | True   |
| `config_profile_contact`     | Institutional config contact information.                                                                                                                                                                                                                                                                                                                                                         | `string` |                                                          |          | True   |
| `config_profile_url`         | Institutional config URL link.                                                                                                                                                                                                                                                                                                                                                                    | `string` |                                                          |          | True   |

## Max job request options

Set the top limit for requested resources for any single job.

| Parameter    | Description                                                                                                                                                                                                                                                                 | Type      | Default | Required | Hidden |
| ------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------- | ------- | -------- | ------ |
| `max_cpus`   | Maximum number of CPUs that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the CPU requirement for each process. Should be an integer e.g. `--max_cpus 1`</small></details>                                      | `integer` | 16      |          | True   |
| `max_memory` | Maximum amount of memory that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the memory requirement for each process. Should be a string in the format integer-unit e.g. `--max_memory '8.GB'`</small></details> | `string`  | 128.GB  |          | True   |
| `max_time`   | Maximum amount of time that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the time requirement for each process. Should be a string in the format integer-unit e.g. `--max_time '2.h'`</small></details>        | `string`  | 240.h   |          | True   |

## Generic options

Less common options for the pipeline, typically set in a config file.

| Parameter                | Description                                                                                                                                                                                                                                                                                                                                                                                                  | Type      | Default                        | Required | Hidden |
| ------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | --------- | ------------------------------ | -------- | ------ |
| `help`                   | Display help text.                                                                                                                                                                                                                                                                                                                                                                                           | `boolean` |                                |          | True   |
| `publish_dir_mode`       | Method used to save pipeline results to output directory. <details><summary>Help</summary><small>The Nextflow `publishDir` option specifies which intermediate files should be saved to the output directory. This option tells the pipeline what method should be used to move these files. See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.</small></details> | `string`  | copy                           |          | True   |
| `email_on_fail`          | Email address for completion summary, only when pipeline fails. <details><summary>Help</summary><small>An email address to send a summary email to when the pipeline is completed - ONLY sent if the pipeline does not exit successfully.</small></details>                                                                                                                                                  | `string`  |                                |          | True   |
| `plaintext_email`        | Send plain-text email instead of HTML.                                                                                                                                                                                                                                                                                                                                                                       | `boolean` |                                |          | True   |
| `max_multiqc_email_size` | File size limit when attaching MultiQC reports to summary emails.                                                                                                                                                                                                                                                                                                                                            | `string`  | 25.MB                          |          | True   |
| `monochrome_logs`        | Do not use coloured log outputs.                                                                                                                                                                                                                                                                                                                                                                             | `boolean` |                                |          | True   |
| `multiqc_config`         | Custom config file to supply to MultiQC.                                                                                                                                                                                                                                                                                                                                                                     | `string`  |                                |          | True   |
| `tracedir`               | Directory to keep pipeline Nextflow logs and reports.                                                                                                                                                                                                                                                                                                                                                        | `string`  | ${params.outdir}/pipeline_info |          | True   |
| `validate_params`        | Boolean whether to validate parameters against the schema at runtime                                                                                                                                                                                                                                                                                                                                         | `boolean` | True                           |          | True   |
| `show_hidden_params`     | Show all params when using `--help` <details><summary>Help</summary><small>By default, parameters set as _hidden_ in the schema are not shown on the command line when a user runs with `--help`. Specifying this option will tell the pipeline to show all parameters.</small></details>                                                                                                                    | `boolean` |                                |          | True   |
| `enable_conda`           | Run this workflow with Conda. You can also use '-profile conda' instead of providing this parameter.                                                                                                                                                                                                                                                                                                         | `boolean` |                                |          | True   |
