## Local, unversioned outputs. Everything here is ignored by Git (see .gitignore).

# outputs/

This directory is the default location for all generated outputs from the analysis scripts in this repository.

In general:
- **Plots, tables, logs, and summary files** are written to `outputs/`
- **Large intermediate artefacts that are reused downstream** (e.g. SCENIC input looms, motif databases) are written to `resources/`

## Subdirectories

### `01_loom_prep/`
Outputs from `scripts/01_loom_prep/`:
- QC plots (e.g. gene/cell filtering summaries)
- Diagnostic figures
- `sessionInfo` files

The SCENIC **input loom** itself is written to `resources/pyscenic_resources/`, not here.

### `02_scenic_multirun/`
Outputs from the SCENIC multirun (Nextflow) pipeline:
- Per-run result folders
- Aggregated SCENIC looms
- Logs and pipeline reports

These outputs are later used by the downstream analysis in `scripts/03_downstream/`.

### `03_downstream/`
Outputs from downstream integration and analysis:
- UMAPs (RNA and SCENIC)
- RSS plots
- Heatmaps
- ORA (GO / Reactome / Hallmark) results
- Exported regulon gene lists
- `sessionInfo` files for reproducibility

## Notes
- This directory is intentionally tracked with a `.gitkeep` but **most contents are ignored by git**.
- Re-running scripts will overwrite files unless timestamped outputs are enabled.
