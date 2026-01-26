# Config File

**Purpose:** Nextflow config for SCENIC multirun (VSN pipelines).

## Primary config
`nf_CPUopt_scenic_multirun_25runs_minGene5.config`, the config used in our analysis

- `numRuns = 25`
- `min_genes_regulon = 5`
- `min_regulon_gene_occurrence = 5`
- `rank_threshold = 5000`, `nes_threshold = 3.0`, `auc_threshold = 0.05`, `max_similarity_fdr = 0.001`, `min_genes = 20`
- `container = 'aertslab/pyscenic_scanpy:0.12.0_1.9.1'`
- `skipReports = true`
- *(tracks optional)* keep `tracksDb` / `tracksAnnotation` **commented** unless using ENCODE tracks

**Resource paths:** point to your `${RESOURCES_DIR}` (default in this repo: `resources/pyscenic_resources/`).  
Motif files are **cisTarget v10 (hg38, ±10 kb around TSS)** and the matching **motif→TF** table.

## Use
```bash
nextflow -C config/nf_CPUopt_scenic_multirun_25runs_minGene5.config \
  run HowieJM/vsn-pipelines -entry scenic -r master
# (add -resume to continue a partial run)
```

## For details, see -> scripts/02_scenic_multirun/README.md
