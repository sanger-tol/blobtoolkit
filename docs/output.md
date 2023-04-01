# Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

The directories comply with Tree of Life's canonical directory structure.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Coverage Analysis Files](#coverage-analysis-files) - Output files from Mosdepth and Count_Busco_Genes
- [BUSCO and diamond_blastp Analysis Files](#blast-analysis-files) - Output files from Diamond_BlastP
- [Blobtools Files](#blobtools-Files) - Output files from blobtools subworkflow
- [Pipeline Information](#pipeline-information) - Report metrics generated during the workflow execution

### Coverage Analysis Files

Here are the files you can expect in the `blobtoolkit/` sub-directory.

```text
blobtoolkit
  ├── mMelMel1_T1.mosdepth.global.dist.txt
  ├── mMelMel1_T1.mosdepth.region.dist.txt
  ├── mMelMel1_T1.mosdepth.summary.txt
  ├── mMelMel1_T1.per-base.bed.gz
  ├── mMelMel1_T1.per-base.bed.gz.csi
  ├── mMelMel1_T1.regions.bed.gz
  ├── mMelMel1_T1.regions.bed.gz.csi
  ├── GCA_922984935.2.subset_busco_genes_count.tsv
  ├── GCA_922984935.2.subset_coverage.tsv
  |── GCA_922984935.2.subset_window_stats.0.1.tsv
  |── GCA_922984935.2.subset_window_stats.0.01.tsv
  |── GCA_922984935.2.subset_window_stats.100000.tsv
  |── GCA_922984935.2.subset_window_stats.1000000.tsv
```

- `*.mosdepth.global.dist.txt` and `*.mosdepth.region.dist.txt` : Text files with global and region cumulative coverage distribution respectively
- `*.summary.txt`: Text file with summary mean depths per chromosome and regions
- `*.per-base.bed.gz` and `*.per-base.bed.gz.csi` : BED file with per-base coverage and it's index file
- `*.regions.bed.gz` and `*.regions.bed.gz.csi` : BED file with per-region coverage and it's index file
- `*_busco_genes_count.tsv` : TSV file counting BUSCOs for different lineages across BED file
- `*_coverage.tsv` : TSV file combining the BED file with per-region coverage from mosdepth and the busco_count_genes TSV file
- `*_window_stats.0.1.tsv` : Aggregate 1kb windows into a window of 10% of the contig length
- `*_window_stats.0.01.tsv` : Aggregate 1kb windows into a window of 1% of the contig length
- `*_window_stats.100000.tsv` : Aggregate 1kb windows into 100kb windows
- `*_window_stats.1000000.tsv` : Aggregate 1kb windows into 1Mb windows

### BUSCO and diamond_blastp Analysis Files

- `blobtoolkit/busco_diamond/`

```text
busco_diamond/
├── GCA_922984935_2-archaea_odb10-busco
├── GCA_922984935_2-archaea_odb10-busco.batch_summary.txt
├── GCA_922984935_2-bacteria_odb10-busco
├── GCA_922984935_2-bacteria_odb10-busco.batch_summary.txt
├── GCA_922984935_2_busco_genes.fasta
├── GCA_922984935_2-carnivora_odb10-busco
├── GCA_922984935_2-carnivora_odb10-busco.batch_summary.txt
├── GCA_922984935_2-eukaryota_odb10-busco
├── GCA_922984935_2-eukaryota_odb10-busco.batch_summary.txt
├── GCA_922984935_2-eutheria_odb10-busco
├── GCA_922984935_2-eutheria_odb10-busco.batch_summary.txt
├── GCA_922984935_2-laurasiatheria_odb10-busco
├── GCA_922984935_2-laurasiatheria_odb10-busco.batch_summary.txt
├── GCA_922984935_2-mammalia_odb10-busco
├── GCA_922984935_2-mammalia_odb10-busco.batch_summary.txt
├── GCA_922984935_2-metazoa_odb10-busco
├── GCA_922984935_2-metazoa_odb10-busco.batch_summary.txt
├── GCA_922984935_2-tetrapoda_odb10-busco
├── GCA_922984935_2-tetrapoda_odb10-busco.batch_summary.txt
├── GCA_922984935_2.tsv
├── GCA_922984935_2.txt
├── GCA_922984935_2-vertebrata_odb10-busco
├── GCA_922984935_2-vertebrata_odb10-busco.batch_summary.txt
├── short_summary.specific.archaea_odb10.GCA_922984935.2.subset.fasta.json
├── short_summary.specific.archaea_odb10.GCA_922984935.2.subset.fasta.txt
├── short_summary.specific.bacteria_odb10.GCA_922984935.2.subset.fasta.json
├── short_summary.specific.bacteria_odb10.GCA_922984935.2.subset.fasta.txt
├── short_summary.specific.carnivora_odb10.GCA_922984935.2.subset.fasta.json
├── short_summary.specific.carnivora_odb10.GCA_922984935.2.subset.fasta.txt
├── short_summary.specific.eukaryota_odb10.GCA_922984935.2.subset.fasta.json
├── short_summary.specific.eukaryota_odb10.GCA_922984935.2.subset.fasta.txt
├── short_summary.specific.eutheria_odb10.GCA_922984935.2.subset.fasta.json
├── short_summary.specific.eutheria_odb10.GCA_922984935.2.subset.fasta.txt
├── short_summary.specific.laurasiatheria_odb10.GCA_922984935.2.subset.fasta.json
├── short_summary.specific.laurasiatheria_odb10.GCA_922984935.2.subset.fasta.txt
├── short_summary.specific.mammalia_odb10.GCA_922984935.2.subset.fasta.json
├── short_summary.specific.mammalia_odb10.GCA_922984935.2.subset.fasta.txt
├── short_summary.specific.metazoa_odb10.GCA_922984935.2.subset.fasta.json
├── short_summary.specific.metazoa_odb10.GCA_922984935.2.subset.fasta.txt
├── short_summary.specific.tetrapoda_odb10.GCA_922984935.2.subset.fasta.json
├── short_summary.specific.tetrapoda_odb10.GCA_922984935.2.subset.fasta.txt
├── short_summary.specific.vertebrata_odb10.GCA_922984935.2.subset.fasta.json
└── short_summary.specific.vertebrata_odb10.GCA_922984935.2.subset.fasta.txt
```

- `GCA_922984935_2.txt`: Text file containing hits in tabular BLAST format.
- `GCA_922984935_2.tsv`: Tab-delimited file containing GOAT taxon search results.
- `*_odb10*`: BUSCO results for each lineage in GOAT taxon search results.

### Blobtools Files

Output config file from blobltools subworkflow is exported to the `blobltoolkit/ACCESSION` folder:

```
GCA_922984935.2
└── config.yaml
```

Output blobdir files are exported to `blobdir/meta` folder (where meta is either the accession or ToLID):

```
GCA_922984935_2/
├── archaea_odb10_count.json
├── archaea_odb10_count_windows_1000000.json
├── archaea_odb10_count_windows_100000.json
├── bacteria_odb10_count.json
├── bacteria_odb10_count_windows_1000000.json
├── bacteria_odb10_count_windows_100000.json
├── buscogenes_class_cindex.json
├── buscogenes_class.json
├── buscogenes_class_positions.json
├── buscogenes_class_score.json
├── buscogenes_class_windows.json
├── buscogenes_family_cindex.json
├── buscogenes_family.json
├── buscogenes_family_positions.json
├── buscogenes_family_score.json
├── buscogenes_family_windows.json
├── buscogenes_genus_cindex.json
├── buscogenes_genus.json
├── buscogenes_genus_positions.json
├── buscogenes_genus_score.json
├── buscogenes_genus_windows.json
├── buscogenes_kingdom_cindex.json
├── buscogenes_kingdom.json
├── buscogenes_kingdom_positions.json
├── buscogenes_kingdom_score.json
├── buscogenes_kingdom_windows.json
├── buscogenes_order_cindex.json
├── buscogenes_order.json
├── buscogenes_order_positions.json
├── buscogenes_order_score.json
├── buscogenes_order_windows.json
├── buscogenes_phylum_cindex.json
├── buscogenes_phylum.json
├── buscogenes_phylum_positions.json
├── buscogenes_phylum_score.json
├── buscogenes_phylum_windows.json
├── buscogenes_positions.json
├── buscogenes_species_cindex.json
├── buscogenes_species.json
├── buscogenes_species_positions.json
├── buscogenes_species_score.json
├── buscogenes_species_windows.json
├── buscogenes_superkingdom_cindex.json
├── buscogenes_superkingdom.json
├── buscogenes_superkingdom_positions.json
├── buscogenes_superkingdom_score.json
├── buscogenes_superkingdom_windows.json
├── carnivora_odb10_busco.json
├── carnivora_odb10_count.json
├── carnivora_odb10_count_windows_1000000.json
├── carnivora_odb10_count_windows_100000.json
├── eukaryota_odb10_count.json
├── eukaryota_odb10_count_windows_1000000.json
├── eukaryota_odb10_count_windows_100000.json
├── eutheria_odb10_count.json
├── eutheria_odb10_count_windows_1000000.json
├── eutheria_odb10_count_windows_100000.json
├── GCA_922984935_2_cov.json
├── GCA_922984935_2_cov_stats.json
├── GCA_922984935_2_cov_windows_1000000.json
├── GCA_922984935_2_cov_windows_100000.json
├── identifiers.json
├── laurasiatheria_odb10_count.json
├── laurasiatheria_odb10_count_windows_1000000.json
├── laurasiatheria_odb10_count_windows_100000.json
├── length.json
├── length_windows_1000000.json
├── length_windows_100000.json
├── mammalia_odb10_count.json
├── mammalia_odb10_count_windows_1000000.json
├── mammalia_odb10_count_windows_100000.json
├── meta.json
├── metazoa_odb10_count.json
├── metazoa_odb10_count_windows_1000000.json
├── metazoa_odb10_count_windows_100000.json
├── position.json
├── position_windows_1000000.json
├── position_windows_100000.json
├── proportion.json
├── proportion_windows_1000000.json
├── proportion_windows_100000.json
├── tetrapoda_odb10_count.json
├── tetrapoda_odb10_count_windows_1000000.json
├── tetrapoda_odb10_count_windows_100000.json
├── vertebrata_odb10_count.json
├── vertebrata_odb10_count_windows_1000000.json
└── vertebrata_odb10_count_windows_100000.json
```

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
