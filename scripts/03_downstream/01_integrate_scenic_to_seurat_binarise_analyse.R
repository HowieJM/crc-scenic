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
                 
# Ensure baseline default assay is RNA
if (!"RNA" %in% names(so@assays)) {
  stop("Assay 'RNA' not found in Seurat object. Available: ",
       paste(names(so@assays), collapse = ", "))
}
DefaultAssay(so) <- "RNA"
                 
## Optional: external TAS annotation CSV (Magdalena’s updated labels, if not already added)
ANNO_CSV <- Sys.getenv("TAS_ANNOT_CSV", unset = "data/Stroma_Annotation_Seurat.csv")
anno <- read.csv(ANNO_CSV, header = TRUE)
str(anno) # these annotations applied to the entire pre-SCENIC filtered Seurat object

# Retain only cells present in the pre-SCENIC filtered Seurat object
seurat_cells  <- so@meta.data$cell.ID
filtered_anno <- anno[anno$cell.ID %in% seurat_cells, ]
cat("Number of matching cells:", nrow(filtered_anno), "\n")  # expect 24044
             
# Add updated annotations to Seurat metadata
rownames(filtered_anno) <- filtered_anno$cell.ID
filtered_anno <- filtered_anno[, -1]  # drop cell.ID column
colnames(filtered_anno) <- c("seurat_clusters_New", "Subcluster_New")  # so we don't lose track
so <- AddMetaData(so, metadata = filtered_anno)
head(so@meta.data[, c("seurat_clusters_New", "Subcluster_New")])

# Final safety check: same cell order
stopifnot(isTRUE(all.equal(rownames(filtered_anno), so@meta.data$cell.ID)))
str(so) # updated TAS annotations added to the pre-SCENIC filtered Seurat object

##

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~          
## Integrate SCENIC outputs into Seurat (add motif-based regulon information)

# Open loom using SCopeLoomR
loom <- open_loom(SCENIC_LOOM, mode = "r")

# Extract regulon incidence matrix (MotifRegulons in row_attrs; which genes per regulon), and export for record
regulonsMat <- SCopeLoomR::get_regulons(loom, column.attr.name = "MotifRegulons")
write.csv(regulonsMat, file = file.path(out_folder, "pySCENIC_regulons_incidence_matrix.csv"), 
         row.names = TRUE)
message("Wrote regulon incidence: ", file.path(out_folder, "pySCENIC_regulons_incidence_matrix.csv"))

# Convert incidence to gene list, and export for record -> list of genes per regulon 
regulons <- SCENIC::regulonsToGeneLists(regulonsMat)
saveRDS(regulons, file = file.path(out_folder, "pySCENIC_regulons_list.rds"))
sink(file.path(out_folder, "pySCENIC_regulons_list.txt"))
for (regulon in names(regulons)) {
  cat("Regulon:", regulon, "\nGenes:", paste(regulons[[regulon]], collapse = ", "), "\n\n")
}
sink() # save list
message("Wrote regulon lists (RDS + TXT) to: ", out_folder)

# Extract regulon activity (AUC) per cell
regulonAUC <- SCopeLoomR::get_regulons_AUC(loom, column.attr.name = "MotifRegulonsAUC",
                                           rows = "regulons", columns = "cells")
AUCmat <- AUCell::getAUC(regulonAUC) # take AUC activity matrix in plain numeric form 

# Add continuous AUC assay to Seurat
so[['pyscenicAUC']] <- CreateAssayObject(data = AUCmat)  # regulon activity is now added
stopifnot(identical(rownames(AUCmat), rownames(so[["pyscenicAUC"]]))) # quick check
str(so)
                 
# Extract thresholds for binarisation: # Note -> the list includes candidate, non-aggregated regulons
regulonAucThresholds <- SCopeLoomR::get_regulon_thresholds(loom, only.selected = TRUE)                
# Also note: in this loom the vector is inverted (names = thresholds, values = regulon IDs); so ->

# Re-build as a named numeric vector: names = regulon, values = threshold
thr_vec <- setNames(
  as.numeric(names(regulonAucThresholds)),          # thresholds (numeric)
  as.character(unlist(regulonAucThresholds))        # regulon IDs (character)
)

# Export CSV in the natural column order
thresholds_df <- data.frame(
  Regulon   = names(thr_vec),
  Threshold = as.numeric(thr_vec),
  row.names = NULL
)
write.csv(thresholds_df,
          file = file.path(out_folder, "regulonAUC_thresholds.csv"),
          row.names = FALSE)
                 

# Binarise the Activity Matrix:

# use regulon activity matrix
str(AUCmat)  # row = regulons, col = cells

# use regulon activity thresholds, in required format (names = regulon, values = cutoff)
common_regulons <- intersect(rownames(AUCmat), names(thr_vec))  # expect aggregated set (e.g. 288)
AUCbin <- 1L * sweep(AUCmat[common_regulons, , drop = FALSE], 1, thr_vec[common_regulons], FUN = ">")  # binarise by row (aka by regulon)

# add binarised assay
so[["pyscenicAUC_bin"]] <- CreateAssayObject(data = AUCbin)
dim(so[["pyscenicAUC_bin"]])  # expect 288 24044, proceed

# We now have a continuous measure of regulon activity, and a binary "active" or "inactive" state per cell              
          

# Add SCENIC embeddings (UMAP/t-SNE computed on AUC in the loom)
embeddings <- get_embeddings(loom)  # list 
stopifnot("SCENIC AUC UMAP" %in% names(embeddings), "SCENIC AUC t-SNE" %in% names(embeddings))

umap_matrix <- embeddings[["SCENIC AUC UMAP"]]
tsne_matrix <- embeddings[["SCENIC AUC t-SNE"]]

# Ensure cell alignment (rows) and name dims for Seurat
if (!identical(rownames(umap_matrix), colnames(so))) {
  # If needed, align by cell name; error if mismatch
  stopifnot(all(rownames(umap_matrix) %in% colnames(so)))
  umap_matrix <- umap_matrix[colnames(so), , drop = FALSE]
}
colnames(umap_matrix) <- c("UMAP_1", "UMAP_2")

if (!identical(rownames(tsne_matrix), colnames(so))) {
  stopifnot(all(rownames(tsne_matrix) %in% colnames(so)))
  tsne_matrix <- tsne_matrix[colnames(so), , drop = FALSE]
}
colnames(tsne_matrix) <- c("tSNE_1", "tSNE_2")

so[["SCENIC_UMAP"]]  <- CreateDimReducObject(embeddings = umap_matrix, key = "SCENICUMAP_",  assay = "pyscenicAUC")
so[["SCENIC_tSNE"]]  <- CreateDimReducObject(embeddings = tsne_matrix, key = "SCENICtSNE_",   assay = "pyscenicAUC")

# Re-label lsTAS -> to match paper
so$Subcluster_New[so$Subcluster_New == "lsTAS"] <- "apTAS"

# Levels
so$Subcluster_New <- factor(so$Subcluster_New,
                            levels = c("apTAS",
                                       "pTAS",
                                       "iTAS1", "iTAS2", 
                                       "mscTAS", 
                                       "myTAS1", "myTAS2", "myTAS3"))

# Save integrated Seurat -> uncomment to save, can reinitiate downstream analysis from here
# saveRDS(so, file = file.path(out_folder,
#   "DATA_Joanito_Stroma-only_filtered_for_SCENIC_input_withBinarizedRegulonActivity_and_Regulon_UMAP_tSNE.rds"))

# Close loom
close_loom(loom)

                 
#

                 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Analyse Data -> ID Cell Type Specific Regulons, etc 
                 

# 1 - Initial Overview: Cell Types on RNA and Regulon Activity UMAPs
u1 <- DimPlot(so, reduction = "umap",        group.by = "Subcluster_New") + ggtitle("RNA-based UMAP")
u2 <- DimPlot(so, reduction = "SCENIC_UMAP", group.by = "Subcluster_New") + ggtitle("SCENIC AUC-based UMAP")
u3 <- DimPlot(so, reduction = "SCENIC_tSNE", group.by = "Subcluster_New") + ggtitle("SCENIC AUC-based t-SNE")

out_F <- file.path(plot_folder, "1-UMAP-CellTypeProjections_RNA_UMAP__SCENIC_UMAP__SCENIC_tSNE.pdf")
pdf(out_F, width = 18, height = 6)
plot(u1 + u2 + u3 + plot_layout(ncol = 3, widths = c(1,1,1)))
dev.off()

# 2 - Identify TAS-specific regulons (continuous AUC) — FindAllMarkers gives *enriched* (top) regulons per TAS,
#    not necessarily unique; we assess specificity later via RSS.                 
DefaultAssay(so) <- "pyscenicAUC"
Idents(so) <- so$Subcluster_New

regulon_markers <- FindAllMarkers(
  object = so,
  assay  = "pyscenicAUC",
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.10
)
write_csv(regulon_markers, file.path(plot_folder, "2-RegulonResults_TopPerCellType_FindMarkers_pct_0_1__logfc_0_1.csv"))

DefaultAssay(so) <- "RNA"   
                 
# Plot top single regulon activity marker per TAS 
top1_regulon_markers <- regulon_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 1, with_ties = FALSE) %>%
  ungroup()
print(top1_regulon_markers, n = 8)
write_csv(top1_regulon_markers, file.path(plot_folder, "2A-Top1_regulon_per_TAS.csv"))

# Optional: quick feature plots of top-1 regulon per TAS on RNA-UMAP and SCENIC UMAP
DefaultAssay(so) <- "pyscenicAUC"

pairs <- setNames(top1_regulon_markers$gene, top1_regulon_markers$cluster)

# RNA UMAP
p_rna <- imap(pairs, ~ FeaturePlot(so, features = .x, reduction = "umap") +
                        labs(title = paste0(.y, ": ", .x)))
pdf(file.path(plot_folder, "2B-Top1_Regulon_Per_TAS_RNAUMAP.pdf"), width = 12, height = 10)
print(wrap_plots(p_rna, ncol = 3))
dev.off()

# SCENIC UMAP
p_scenic <- imap(pairs, ~ FeaturePlot(so, features = .x, reduction = "SCENIC_UMAP") +
                           labs(title = paste0(.y, ": ", .x)))
pdf(file.path(plot_folder, "2C-Top1_Regulon_Per_TAS_SCENICUMAP.pdf"), width = 12, height = 10)
print(wrap_plots(p_scenic, ncol = 3))
dev.off()

DefaultAssay(so) <- "RNA"

# Optional: repeat top-1 regulon plots using binarised activity (on/off)
# DefaultAssay(so) <- "pyscenicAUC_bin"
# p_rna_bin <- imap(pairs, ~ FeaturePlot(so, features = .x, reduction = "umap") +
#                             labs(title = paste0(.y, " (BIN): ", .x)))
# pdf(file.path(plot_folder, "2D-Top1_Regulon_Per_TAS_RNAUMAP_BIN.pdf"), width = 12, height = 10)
# print(wrap_plots(p_rna_bin, ncol = 3))
# dev.off()
# DefaultAssay(so) <- "RNA"
                 
# Visual QC / first-pass interpretation:
# - Top regulons per TAS show expected correspondence with TAS clusters;
# - Patterns are similar on RNA vs SCENIC embeddings (different trait spaces, similar separation);
# - Specificity varies by regulon and TAS (gradients are visible).

#


                 
                 
                 
                 




