# ======================================================================================
# Title: Import SCENIC multirun loom → integrate with Seurat → binarise regulon activity
# Author: James Malcolm Howie
# Repository: https://github.com/HowieJM/crc-scenic
#
# Purpose
#   Read the aggregated SCENIC loom (multirun), add the AUC assays & SCENIC embeddings
#   to the initial TAS Seurat object, and create a binarised regulon-activity matrix.
#
#   Then, run downstream analysis of regulon activity (RSS) and pathway (ORA) analysis
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

# Optional: override via env vars without editing the script
SCENIC_LOOM   <- Sys.getenv(
  "SCENIC_LOOM",
  unset = "outputs/02_scenic_multirun/<RUN_TAG>/out/scenic/<project>__25Runs__YYYY_MM_DD/data/<project>__25Runs__YYYY_MM_DD.SCENIC.loom"
)
TAS_SEURAT_RDS <- Sys.getenv(
  "TAS_SEURAT_RDS",
  unset = "data/TAS_seurat.rds"   # not bundled; set to your TAS Seurat path
)
                 
                 
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


#~~~~~~~~~~~~~~~~~~~~~
## Integrate Data Sets
                 
# Either run the below sections, or skip to the analysis step, and load the pre-build so, with binarised regulon activity matrix


#~~~~~~~~~~~~~~~~~~~~~
## Inspect SCENIC Data - optional quick checks: VSN-pySCENIC MultiRun Loom 

# Read SCENIC out loom
loom_file <- H5File$new(SCENIC_LOOM, mode = "r")
                 
# Minimal non-interactive checks ->

# Loom spec version
loom_spec <- loom_file[["attrs/LOOM_SPEC_VERSION"]][]
message("LOOM_SPEC_VERSION: ", loom_spec)  # expect 3.0.0

# Gene-Cell matrix size
genes <- nrow(loom_file[["matrix"]][] ) 
cells <- ncol(loom_file[["matrix"]][] )
message("Genes × Cells: ", genes, " × ", cells)  # expect 24719 genes, 24,044 cells

# Number of stable aggregated regulons (AUC; ≥80% of runs)    
n_regs <- ncol(loom_file[["col_attrs/MotifRegulonsAUC"]][])
message("Stable Regulons (AUC; ≥80% of runs): ", n_regs)  #expect 288 regulons 


# Additional deeper checks (unhash if desired):
#loom_file$ls(recursive = TRUE)
#loom_file[["layers"]]

#loom_file[["col_attrs"]] #cell metadata
#loom_file[["row_attrs"]] #gene metadata

#loom_file[["col_attrs/CellID"]]                         
#loom_file[["row_attrs/Gene"]]                           
#loom_file[["col_attrs/MotifRegulonsAUC"]]   
#loom_file[["row_attrs/MotifRegulonGeneOccurrences"]]   

#loom_file[["row_attrs/MotifRegulonGeneWeights"]]
#loom_file[["row_attrs/MotifRegulons"]]

#loom_file[["matrix"]] #expression matrix

#ncol(head(loom_file[["col_attrs/MotifRegulonsAUC"]][])) # number of regulons
#ncol(head(loom_file[["row_attrs/MotifRegulonGeneOccurrences"]][])) # number of regulons

                 
# Close the Loom file
loom_file$close_all()



#~~~~~~~~~~~~~~~~~~
## Load Seurat data

# TAS Seurat -> Data provenance (summary, see scripts/01_loom_prep/README.md for full):

# - scRNA-seq: Joanito et al., Nat Genet 2022 (CRC; 10x Genomics; GRCh38)
# - TAS subset: >20 stromal cells per patient, all cohorts (M. Frank re-annotated TAS Seurat)
# - filtered for SCENIC (J. M. Howie, see scripts/01_loom_prep/01_tas_seurat_filtered_to_scenic_loom.R)
     
so <- readRDS(TAS_SEURAT_RDS)
str(so) 
                 
                 

                 
                 
                 




