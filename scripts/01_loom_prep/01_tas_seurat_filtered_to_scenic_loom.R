# =============================================================================================
# Title: TAS-only scRNA Seurat filtering → export loom for SCENIC
# Author: James Malcolm Howie
# Repository: https://github.com/HowieJM/crc-scenic
#
# Provenance:
#   - Final:   9-20240902_PrepFilteredscRNA_for_PyScenic-ModifedExport.R
#
# Purpose: filter scRNA Seurat data and create loom file → export SCENIC-ready loom.
#
# Details:
# Take scRNA TAS-only subset across all CRC cohorts from Joanito et al. (Nat Genet 2022);
# Filter TAS-only Seurat based on SCENIC requirements, then make and export SCENIC-ready loom.
# Filters used: nFeature_RNA 1250–5000; percent.mt ≤ 7.5%; keep genes seen in ≥3 cells.
# =============================================================================================


#~~~~~~~~~~~~~~~~~~~~
# RStudio v 2023.09.0
# R v 4.4.0

#

#~~~~~~~~~~~~~~~~~~~~
# Prepare Environment

out_folder <- file.path("outputs", "20241105_PrepFilteredscRNA_for_PyScenic_to_save_seurat_RDS")
dir.create(out_folder, recursive = TRUE, showWarnings = FALSE) #for plots

resources_dir <- Sys.getenv("RESOURCES_DIR", unset = file.path("resources", "pyscenic_resources"))
dir.create(resources_dir, recursive = TRUE, showWarnings = FALSE) #for loom, to use in SCENIC

set.seed(1414)


#~~~~~~~~~~~~~~
# Load Packages

# Non-standard install
#install.packages("remotes")
#remotes::install_github("aertslab/SCopeLoomR")

# Packages
suppressPackageStartupMessages({
  library(Seurat)      # 4.4.0
  library(tidyverse)
  library(SCopeLoomR)  # loom write/inspect
  library(hdf5r)       # optional: used for H5File() inspection
})


#~~~~~~~~~~
# Load Data 

# Seurat Object - Data -> scRNA, Joanito et al.2022, colorectal, 10X GRCh38 - stromal cells - patients with >20 stromal cells, 5 cohorts
so <- readRDS("DATA_Joanito_Stroma-only.rds") #contact Dietmar if requiring file 

DefaultAssay(so) <- "RNA"   #to extract raw counts



#~~~~~~~
# Filter 

# Factorise Cohort
var_levels <- list()
var_levels[["cohort"]] <- c("CRC-SG1", "CRC-SG2", "KUL3", "KUL5", "SMC")
so@meta.data$cohort <- factor(so@meta.data$dataset, levels = var_levels[["cohort"]])

v1<-VlnPlot(so,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        pt.size = 0,
        group.by = "cohort")
v1

out_F <- paste(out_folder, "/seurat.qc.png", sep = "")
png(out_F, width = 1000, height = 400)
plot(v1)
dev.off()


# Filter low-high expressed genes, and high mitochondrial content cells
nFeature_RNA_min <- 1250   # empirical, features per cell
nFeature_RNA_max <- 5000   # ""
percent.mt_max   <- 7.5    # recommended 5-15%, 10% fits 

# I wanted to be relatively permissive, but cut low quality cells; thresholds cut similar inflection to their tutorial #https://github.com/aertslab/SCENICprotocol/blob/master/notebooks/PBMC10k_SCENIC-protocol-CLI.ipynb

dim(so) #genes, cells
so <- subset(so, nFeature_RNA >= nFeature_RNA_min & nFeature_RNA <= nFeature_RNA_max & percent.mt <= percent.mt_max)    
dim(so)

v2<-VlnPlot(so,
            features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
            pt.size = 0,
            group.by = "cohort")
v2

out_F <- paste(out_folder, "/seurat.qc.filtered.png", sep = "")
png(out_F, width = 1000, height = 400)
plot(v2)
dev.off()


# Filter out genes expressed in very few cells
min_cells_expressed <- 3   # recommended

gene_counts <- Matrix::rowSums(so@assays$RNA@counts > 0)
genes_to_keep <- names(gene_counts)[gene_counts >= min_cells_expressed]

dim(so) #genes, cells
so <- subset(so, features = genes_to_keep)
dim(so) 

v3<-VlnPlot(so,
            features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
            pt.size = 0,
            group.by = "cohort")
v3

out_F <- paste(out_folder, "/seurat.qc.filtered_min_cell_filtered.png", sep = "")
png(out_F, width = 1000, height = 400)
plot(v3)
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Check Data + Deal with NAs etc

# Check metadata integrity
check_metadata_na <- function(seurat_object) {
  # Check for NA values in the metadata
  if (any(is.na(seurat_object@meta.data))) {
    cat("Warning: NA values found in metadata.\n")
    
    # Identify columns with NA values and their counts
    na_columns <- colnames(seurat_object@meta.data)[colSums(is.na(seurat_object@meta.data)) > 0]
    na_counts <- colSums(is.na(seurat_object@meta.data))
    
    # Print columns with NA values and their counts
    cat("Columns with NA values and their counts:\n")
    for (col in na_columns) {
      cat(col, ": ", na_counts[col], " NA values\n", sep = "")
    }
  } else {
    cat("No NA values found in metadata.\n")
  }
}
check_metadata_na(so)

# iCMS and msi are factors with one level, mostly NA, hence we remove these for now [can re-add later from reference csv files] 
cols_to_remove <- c("iCMS", "msi")
so@meta.data <- so@meta.data %>% select(-all_of(cols_to_remove))

# Convert metadata character columns to factors, keeping critical identifiers as characters
so@meta.data <- so@meta.data %>%
  mutate(across(where(is.character) & !c("cell.ID"), as.factor))

# Check dimensions data and metadata match
check_dimensions <- function(seurat_object) {
  # Number of cells and genes
  num_cells <- ncol(seurat_object)
  num_genes <- nrow(seurat_object)
  
  # Number of metadata rows
  num_meta_rows <- nrow(seurat_object@meta.data)
  
  # Print the dimensions
  cat("Number of cells: ", num_cells, "\n")
  cat("Number of genes: ", num_genes, "\n")
  
  # Check if the number of metadata rows matches the number of cells
  if (num_meta_rows == num_cells) {
    cat("Metadata rows match the number of cells.\n")
  } else {
    cat("Warning: Metadata rows (", num_meta_rows, ") do not match the number of cells (", num_cells, ").\n", sep = "")
  }
}
check_dimensions(so)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract Filtered Raw Count Matrix for PyScenic - Retain Metadata

# 1 - Raw, Full Dataset

if (DefaultAssay(so) == "RNA") {
  cat("The default assay is set to RNA, proceed.\n")
} else {
  cat("Caution: The default assay is not set to RNA. It is set to:", DefaultAssay(so), "\n")
}

# Remove integration and dimensional reductions from the Seurat object
sox <- DietSeurat(so, assays ="RNA", dimreducs = FALSE)
print(sox)


# SCopeLoomR: to Produce Loom File   #further notes at: https://github.com/hbc/knowledgebase/blob/master/scrnaseq/pySCENIC.md

# Extract raw count matrix [can also use normalized if integrated retaining all genes]
dgem <- GetAssayData(sox, assay = "RNA", slot = "counts")
dim(dgem)
head(colnames(dgem)) # columns => cell IDs

# Extract cell-level metadata
cell.info <- sox@meta.data

# Extract Gene Metadata (if any)
gene.info <- data.frame(Gene = rownames(dgem))

# Extract default embedding (e.g. UMAP or PCA coordinates)   #if desired; note -> extracting from the non-skinny seurat object    
default.umap <- Embeddings(so, reduction = "umap")
default.umap.name <- "UMAP"

# Create the loom file
file.name <- file.path(
  resources_dir,
  "Joanito_StromaOnly_Filtered_Seurat_Object_All_Cohorts_FullDataSet_SCopeLoomR.loom"
)

project.title <- "Magdalena Stroma Project - Filtered Data, All Cohorts"

build_loom(
  file.name = file.name,
  dgem = dgem,
  title = project.title,
  genome = "human",
  default.embedding = default.umap,
  default.embedding.name = default.umap.name
)


# Inspect Loom:
stopifnot(file.exists(file.name))
loom_file <- H5File$new(file.name, mode = "r")

# List contents
loom_file$ls(recursive = TRUE)

loom_file[["attrs/LOOM_SPEC_VERSION"]][] #3.0.0 -> this is a key global attribute, which MUST be so

loom_file[["layers"]]

loom_file[["col_attrs/CellID"]]
loom_file[["row_attrs/Gene"]]

loom_file[["col_attrs"]] #cell metadata
loom_file[["row_attrs"]] #gene metadata

loom_file[["matrix"]]    #expression matrix

# Close the Loom file
loom_file$close_all()

#

#############################################################
sI <- sessionInfo()
out_file <- paste(out_folder, "/sessionInfo.rds", sep = "")
saveRDS(object = sI, file = out_file)
out_file_txt <- paste(out_folder, "/sessionInfo.txt", sep = "")
writeLines(capture.output(sessionInfo()), con = out_file_txt)
#############################################################
