# 01 — Loom prep

**Purpose:** Build the TAS-only SCENIC input loom from the TAS-only Seurat object, then (optionally) sanity-check it.

## Files in this folder
- `01_tas_seurat_filtered_to_scenic_loom.R` — **Seurat → loom** (to produce SCENIC-ready file)
- `02_check_loom_all.py` — tiny Python checker (version / global attributes)
- `03_jupyter_tunnel.md` — how to open a remote Jupyter via SSH tunnel (to view the loom)

## Inputs
- TAS Seurat RDS (TAS re-annotation from the scRNA analysis in Frank et al. 2026)

## Outputs
- SCENIC input loom → `resources/pyscenic_resources/Joanito_StromaOnly_Filtered_Seurat_Object_All_Cohorts_FullDataSet_SCopeLoomR.loom`
- QC plots + session info → `outputs/01_loom_prep/`

## Quick start
```bash
# build loom (from repo root)
Rscript scripts/01_loom_prep/01_tas_seurat_filtered_to_scenic_loom.R
```

## Optional checks
### print LOOM_SPEC_VERSION
```bash
python scripts/01_loom_prep/02_check_loom_all.py --version resources/pyscenic_resources/Joanito_StromaOnly_Filtered_Seurat_Object_All_Cohorts_FullDataSet_SCopeLoomR.loom
```

### list global attributes
```bash
python scripts/01_loom_prep/02_check_loom_all.py --attrs   resources/pyscenic_resources/Joanito_StromaOnly_Filtered_Seurat_Object_All_Cohorts_FullDataSet_SCopeLoomR.loom
```

