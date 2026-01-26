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
#
# Tested in: RStudio v 2023.09.0  # R v 4.4.0 / 4.4.2, and
#            RStudio v 2024.04.1  # R v 4.4.0
# ======================================================================================



#~~~~~~~~~~~~~~~~~~~~
# Prepare Environment

out_folder  <- file.path("outputs", "03_downstream")
plot_folder <- file.path(out_folder, "plots")
dirs <- c(out_folder, plot_folder)
invisible(lapply(dirs, function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)))

set.seed(1414)


                 
#~~~~~~~~~~~~~~
# Load Packages

# Non-standard install
#install.packages("remotes")
#remotes::install_github("aertslab/SCopeLoomR")
#BiocManager::install("RcisTarget") #in R 4.4.0+, or fails SCENIC
#install.packages("xml2")
#remotes::install_github("aertslab/SCENIC")
#BiocManager::install("AUCell")

# Packages ---
library(Seurat)      # 4.4.0
library(tidyverse)
library(data.table)
library(hdf5r)
library(SCopeLoomR)
library(SCENIC)
library(patchwork)
library(reshape2)
library(dplyr); library(stringr); library(tibble)
library(Matrix)      # for rowMeans on dgCMatrix

                 
# ORA (loaded later when needed)
# library(clusterProfiler); library(org.Hs.eg.db); library(ReactomePA); library(msigdbr)
# library(enrichplot); library(GOSemSim); library(readr)





