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
library(ggtext)

                 
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

                 
# Explore wider range of markers per TAS -> which regulon sets classify each TAS
DefaultAssay(so) <- "pyscenicAUC"

# Top 10 regulon markers per TAS (continuous AUC)
top10_regulon_markers <- regulon_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
  ungroup()

# Features used downstream
feats_top <- unique(top10_regulon_markers$gene)  # unique regulon list
# Note: in current data no TAS has >10 regulons, so this captures all per TAS

# Heatmap (unscaled AUC)
h_unscaled <- DoHeatmap(
  so,
  features = feats_top,
  assay    = "pyscenicAUC",
  slot     = "data",
  raster   = TRUE
) + ggtitle("Top regulon markers by TAS (AUC, unscaled)")

# Scale AUC for these features (row-wise z within regulon)
# Note: ScaleData() modifies the 'pyscenicAUC' assay in 'so' for these features.
so <- ScaleData(
  so,
  assay    = "pyscenicAUC",
  features = feats_top,
  verbose  = FALSE
)

# Heatmap (scaled AUC)
h_scaled <- DoHeatmap(
  so,
  features = feats_top,
  assay    = "pyscenicAUC",
  slot     = "scale.data",
  raster   = TRUE
) +
  ggtitle("Top 10 regulon markers by TAS (scaled AUC)") +
  labs(fill = "Scaled AUC (z)")

# Save PNGs (and optionally PDF)
ggsave(file.path(plot_folder, "7A-Heatmap_TopRegulonMarkers_unscaled.png"),
       h_unscaled, width = 8, height = 10, dpi = 300)
ggsave(file.path(plot_folder, "7B-Heatmap_TopRegulonMarkers_scaled.png"),
       h_scaled,   width = 8, height = 10, dpi = 300)
# ggsave(file.path(plot_folder, "7B-Heatmap_TopRegulonMarkers_scaled.pdf"), h_scaled, width = 8, height = 10)

# Export the table used
write_csv(top10_regulon_markers, file.path(plot_folder, "8-top10_regulon_markers.csv"))

# Interpretation (heatmap of top regulon markers):
# - Clear diagonal: TAS-specific regulon activity is visible, with some overlap across TAS.
# - Counts vary by TAS (not every TAS necessarily has 10 enriched regulons).
# - Specificity differs across regulons/TAS; gradients are common.
# - These “top” regulons (FindAllMarkers) need not be unique; we assess specificity via RSS below.

# Optional: export genes for top-1 regulon per TAS
regList <- regulons[ intersect(names(regulons), top1_regulon_markers$gene) ]
fn <- file.path(plot_folder, "2D-Top1_regulon_per_TAS_genes.txt")
con <- file(fn, "w")
for (rg in names(regList)) {
  writeLines(paste0("Regulon: ", rg, "\nGenes: ", paste(regList[[rg]], collapse = ", "), "\n"), con)
}
close(con)
message("Wrote: ", fn)

                 
#


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## RSS – regulon specificity (continuous AUC)

# The above analysis uses FindAllMarkers (enriched, not unique);
# Now, we use RSS to assess specificity across TAS.

# Check alignment
stopifnot(isTRUE(all.equal(names(so$Subcluster_New), colnames(so[["pyscenicAUC"]]))))

# Calculate RSS (regulon × TAS, continuous AUC)
rss_res <- calcRSS(
  AUC            = as.matrix(GetAssayData(so[["pyscenicAUC"]], slot = "data")),
  cellAnnotation = so$Subcluster_New,
  cellTypes      = levels(so$Subcluster_New)
)
str(rss_res)

                 
# Plot RSS across TAS 

# First, [S.I.] for all regulons × all TAS (no filtering)
rssPlot_all <- plotRSS(
  rss_res,
  zThreshold      = -Inf,   # keep everything
  thr             = -Inf,   # do NOT discard on raw RSS
  cluster_columns = TRUE,
  order_rows      = TRUE
)
ggsave(filename = file.path(plot_folder, "9A-RSS_full_288x8.png"),
       plot     = rssPlot_all$plot, width = 10, height = 32, dpi = 300)

# Second, [MS] for all top RSS regulons with z ≥ 2.5 plus forced myTAS1 regulons
# full-matrix z-score (row-wise within each TAS)
rssNorm_all <- scale(rss_res)                   # 288 × 8
rssNorm_all[rssNorm_all < 0] <- 0               # keep positive part only

# rows to display (z-threshold + force-include myTAS1)
keep_vec  <- c("BACH2-(+)-motif", "ZBED1-(+)-motif", "HMGA2-(+)-motif", "LEF1-(+)-motif")
auto_rows <- rownames(rss_res)[ rowSums(rssNorm_all >= 2.5) > 0 ]
rows_show <- union(keep_vec, auto_rows)

# build dot-heatmap data-frame (raw RSS size; z colour)
rss_sub <- rss_res     [rows_show, , drop = FALSE]
z_sub   <- rssNorm_all [rows_show, , drop = FALSE]

df_raw <- reshape2::melt(rss_sub); colnames(df_raw) <- c("Topic","cellType","RSS")
df_z   <- reshape2::melt(z_sub  ); colnames(df_z)   <- c("Topic","cellType","Z")
rss_df <- merge(df_raw, df_z)

# keep same row order as plotRSS() would give
rowOrder <- rev(
  SCENIC:::.plotRSS_heatmap(z_sub, thr = -Inf,
                            cluster_columns = FALSE,
                            order_rows      = TRUE,
                            verbose = FALSE)@row_names_param$labels)
rss_df$Topic <- factor(rss_df$Topic, levels = rowOrder)

# plot with publication axis orientation: TAS on X, regulons on Y
p_pub <- SCENIC:::dotHeatmap(
  rss_df,
  var.x   = "cellType",   # TAS along x
  var.y   = "Topic",      # regulons along y
  var.size= "RSS",        min.size = 0.5, max.size = 5,
  var.col = "Z",
  col.low = "grey90",
  col.mid = "darkolivegreen3",
  col.high= "darkgreen"
) +
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1))

ggsave(file.path(plot_folder, "9B-RSS_zThres2.5_Top_MyTAS1_forced.png"),
       p_pub, width = 6, height = 8, dpi = 300)


# Third, [optional] (prettier) 9B variant: TAS Y, regulon X
# TAS labels colour blocked; 
# regulon names coloured by their “owner” TAS (max z; myTAS1 forced);

# Set TAS colours
tas_cols <- c(
  "myTAS1"="#0589B6","myTAS2"="#B6D1DA","myTAS3"="#094E95",
  "mscTAS"="#F49A16","iTAS1" ="#8F1B1D","iTAS2" ="#EA5A5D",
  "apTAS" ="#8C3459","pTAS"  ="#A19217"
)

# ensure tas_cols order matches TAS order in the plot
_use <- tas_cols[levels(rss_df$cellType)]
stopifnot(!any(is.na(tas_cols_use)))  # all TAS levels must have colours
                 
# Base dot-heatmap with regulons on X, TAS on Y
p_rev <- SCENIC:::dotHeatmap(
  rss_df,
  var.x    = "Topic",     # regulons on X
  var.y    = "cellType",  # TAS on Y
  var.size = "RSS",       min.size = 0.5, max.size = 5,
  var.col  = "Z",
  col.low  = "grey90",
  col.mid  = "darkolivegreen3",
  col.high = "darkgreen"
)

# --- owner-coloured X labels (regulon names) ---
x_full     <- levels(rss_df$Topic)                          # full regulon names
x_stripped <- sub("-\\(\\+\\)-motif$", "", x_full)          # display-only
keeper_force <- c("BACH2-(+)-motif","ZBED1-(+)-motif","HMGA2-(+)-motif","LEF1-(+)-motif")

stopifnot(all(x_full %in% rownames(rssNorm_all)))           # guard: rows exist in z-matrix
reg2tas_x <- vapply(x_full, function(full){
  if (full %in% keeper_force) return("myTAS1")
  zs <- rssNorm_all[full, ]
  names(zs)[which.max(zs)]
}, character(1))

x_labs_col   <- sprintf("<span style='color:%s'>%s</span>",
                        tas_cols[reg2tas_x], x_stripped)

p_rev_col <- p_rev +
  scale_x_discrete(labels = x_labs_col) +
  theme(
    axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin  = margin(l = 2)
  )

# --- TAS label text & coloured strip on the left ---
ylab_df <- data.frame(
  cellType = factor(levels(rss_df$cellType), levels = levels(rss_df$cellType)),
  x = 1
)

ylab_text <- ggplot(ylab_df, aes(x = 1, y = cellType, label = cellType)) +
  geom_text(hjust = 0.5, size = 4.5) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(r = 2))

ylab_strip <- ggplot(ylab_df, aes(x = 1, y = cellType, fill = cellType)) +
  geom_tile(width = 0.5, height = 0.9) +
  scale_fill_manual(values = tas_cols_use, guide = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(r = 0))

# Compose:  text | strip | dot-plot
p_final <- ylab_text + ylab_strip + p_rev_col +
  patchwork::plot_layout(widths = c(0.12, 0.02, 1))

ggsave(file.path(plot_folder, "9B-RSS_zThres2.5_Top_colouredByZowner.png"),
       p_final, width = 10, height = 5, dpi = 300)



                 # below under dev







                 

# 10 — export top regulons’ gene sets (as defined by rows_show)
regs_subset <- regulons[ intersect(names(regulons), rows_show) ]
sub_txt <- file.path(plot_folder, "10-pySCENIC_regulons_list_RSSz_2.5_myTAS1_inc_subset.txt")
sink(sub_txt)
for (rg in names(regs_subset)) {
  cat("Regulon:", rg, "\nGenes:", paste(regs_subset[[rg]], collapse = ", "), "\n\n")
}
sink(); message("Saved: ", sub_txt)

# (NGFR aside removed)

# (Optional) 9C — AUC-coloured dot-heatmap (size=RSS, colour=AUC), if you like to keep it
# This visualises average AUC per TAS for rows_show in addition to RSS size.
# DefaultAssay(so) <- "pyscenicAUC"
# auc_mat    <- GetAssayData(so[["pyscenicAUC"]], slot = "data")
# cell_order <- colnames(rss_res)
# avgAUC <- sapply(cell_order, function(ct) rowMeans(auc_mat[rows_show, so$Subcluster_New == ct, drop = FALSE]))
# colnames(avgAUC) <- cell_order
# df_auc <- reshape2::melt(avgAUC); colnames(df_auc) <- c("Topic","cellType","AUC")
# df_rss <- reshape2::melt(rss_sub); colnames(df_rss) <- c("Topic","cellType","RSS")
# rss_df_auc <- merge(df_rss, df_auc)
# p_auc <- SCENIC:::dotHeatmap(rss_df_auc, var.x="cellType", var.y="Topic",
#                              var.size="RSS", var.col="AUC",
#                              col.low="grey90", col.mid="lightcoral", col.high="red4")
# ggsave(file.path(plot_folder, "9C-RSS_size_rawAUC_zThres2.5_Top.png"),
#        p_auc, width = 10, height = 5, dpi = 300)

# 12A/12B — heatmaps for rows_show regulons (raw and scaled AUC)
rss_feats <- intersect(rows_show, rownames(so[["pyscenicAUC"]]))  # expect ~35
message(length(rss_feats), " regulons will be plotted.")

# raw AUC
h_raw <- DoHeatmap(so, features = rss_feats, slot = "data", raster = TRUE,
                   group.colors = tas_cols) + ggtitle("RSS-selected regulons – raw AUC")
ggsave(file.path(plot_folder, "12A-Heatmap_RSSselected_rawAUC.png"), h_raw, width = 7, height = 9, dpi = 300)

# scaled AUC
so <- ScaleData(so, assay = "pyscenicAUC", features = rss_feats, verbose = FALSE)
row_labs <- sub("-\\(\\+\\)-motif$", "", rss_feats); names(row_labs) <- rss_feats
h_scaled <- DoHeatmap(so, features = rss_feats, slot = "scale.data", raster = TRUE,
                      group.colors = tas_cols) +
  labs(fill = "Scaled AUC (z)") +
  scale_y_discrete(labels = row_labs) +
  theme(axis.text.y  = element_text(size = 10),
        legend.text  = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.position = "right") +
  guides(colour = "none",
         fill   = guide_colourbar(barwidth = 0.5, barheight = 7))
ggsave(file.path(plot_folder, "12B-Heatmap_RSSselected_scaledAUC.png"),
       h_scaled, width = 12, height = 8, dpi = 300)

# 13/14/14B — lead regulon UMAPs are already above (keep those as-is).
                 
                 
                 
                 
                 




