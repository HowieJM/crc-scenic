[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18471338.svg)](https://doi.org/10.5281/zenodo.18471338)
# crc-scenic

SCENIC-based regulatory network analysis of tumour-associated stroma (TAS) in colorectal cancer 

For - Frank et al 2026      ->      Full Paper at: https://doi.org/10.1016/j.canlet.2025.218228

## Citation

If you use this repository, please cite the Zenodo archive below, *in addition to* the associated Cancer Letters paper:

https://doi.org/10.5281/zenodo.18471338

Formatted citations (BibTeX, APA, etc.) are available via the Zenodo record and GitHub’s “Cite this repository” feature.

## Overview

This repository implements a three-stage workflow:

1. **Prepare SCENIC input loom**  
   (`scripts/01_loom_prep/`)  
   - Filter and QC a TAS Seurat object  
   - Generate a SCENIC-compatible loom  
   - Write reusable inputs to `resources/`

2. **Run SCENIC (multirun)**  
   (`scripts/02_scenic_multirun/`)  
   - Launch SCENIC via a Nextflow pipeline  
   - Aggregate results across runs  
   - Produce stable regulons and activity scores

3. **Downstream integration and analysis**  
   (`scripts/03_downstream/`)  
   - Integrate SCENIC results into Seurat  
   - Visualise regulon activity  
   - Identify TAS-specific regulons (RSS)  
   - Perform pathway enrichment (ORA)

## Directory structure (simplified)
```bash
config/      # Nextflow configuration for SCENIC runs
resources/   # Required inputs (SCENIC looms, motif databases)
outputs/     # Local, unversioned outputs (figures, tables, logs)
scripts/
  ├─ 01_loom_prep/
  ├─ 02_scenic_multirun/
  └─ 03_downstream/
```

## Input data
This pipeline starts from a **pre-filtered, TAS-only Seurat object** derived from colorectal cancer scRNA-seq data.
Detailed data provenance (raw data source, TAS re-annotation, and SCENIC pre-filtering) is documented in `resources/README.md`.

The SCENIC multirun stage is executed via a Nextflow pipeline (see `scripts/02_scenic_multirun/`), based on a fork of the VSN pipelines:
https://github.com/HowieJM/vsn-pipelines

# JMH, 03/02/2026
