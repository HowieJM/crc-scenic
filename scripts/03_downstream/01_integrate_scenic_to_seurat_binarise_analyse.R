# ======================================================================================
# Title: Import SCENIC multirun loom → integrate with Seurat → binarise regulon activity
# Author: James Malcolm Howie
# Repository: https://github.com/HowieJM/crc-scenic
#
# Purpose
#   Read the aggregated SCENIC loom (multirun), add the AUC assays & SCENIC embeddings
#   to the initial TAS Seurat object, and create a binarised regulon-activity matrix.
#
# Inputs
#   - Aggregated SCENIC loom (step 02 output; 
#                   e.g. outputs/02_scenic_multirun/<RUN_TAG>/out/scenic/...SCENIC.loom)
#   - TAS Seurat RDS (from TAS re-annotation; 
#                   e.g. data/TAS_seurat.rds; contact Dietmar if needed)
#
# Outputs
#   - Seurat with 'pyscenicAUC' and 'pyscenicAUC_bin' assays
#   - 'SCENIC_UMAP' and 'SCENIC_tSNE' DimReducs
#   - CSVs: regulon incidence, thresholds, lists
#   - Plots/logs under outputs/03_downstream/
# ======================================================================================

