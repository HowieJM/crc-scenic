# scripts/03_downstream/

Downstream integration and analysis of SCENIC results.

## Purpose

This step:
1. Loads the **aggregated SCENIC multirun loom** produced by `scripts/02_scenic_multirun/`
2. Integrates regulon activity (AUC + binarised states) into a TAS Seurat object
3. Adds SCENIC embeddings (UMAP / t-SNE)
4. Identifies TAS-specific regulons using:
   - differential activity (FindAllMarkers)
   - regulon specificity scores (RSS)
5. Performs pathway enrichment (ORA) on selected regulon gene sets


## Notes on reproducibility

The downstream script provided here reflects the analysis used for the figures and tables reported in Frank et al. (2026), with non-essential exploratory and project-specific code removed for clarity. The workflow has been simplified to focus on the analyses presented in the manuscript. While the logic and parameters are unchanged, minor differences in plotting or intermediate outputs may occur depending on environment, package versions, or input formatting retained or used.

## Main script

### `01_integrate_binarize_analyse.R`

**Inputs**
- Aggregated SCENIC loom  
  (from `outputs/02_scenic_multirun/<RUN_TAG>/out/scenic/...SCENIC.loom`)
- TAS Seurat object (RDS)
- TAS annotation CSV (if not already present in the Seurat object)

**Outputs**
Written to `outputs/03_downstream/`:
- Integrated Seurat object (optional RDS)
- UMAPs (RNA-based and SCENIC-based)
- RSS heatmaps and dot plots
- ORA results (GO:BP, Reactome, Hallmark)
- Exported regulon gene lists
- `sessionInfo` files for reproducibility

## Notes
- This script assumes SCENIC has already been run successfully.
- Heavy intermediate artefacts are **not** written here.
- The ORA step is run with default background databases by design (see comments in script).
