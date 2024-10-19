# Install packages if not already installed
if (!require("Seurat")) install.packages("Seurat")
if (!require("scCATCH")) install.packages("scCATCH")
if (!require("biomaRt")) install.packages("biomaRt")
if (!require("dplyr")) install.packages("dplyr")
if (!require("stringr")) install.packages("stringr")

# Load libraries
library(Seurat)
library(scCATCH)
library(biomaRt)
library(dplyr)
library(stringr)

check_files_exist <- function(dir) {
  required_files <- c("matrix.mtx", "features.tsv", "barcodes.tsv")
  existing_files <- list.files(path = dir)
  
  missing_files <- setdiff(required_files, existing_files)
  
  if (length(missing_files) == 0) {
    message("All required files are present.")
  } else {
    stop(message("The following files are missing: ", paste(missing_files, collapse = ", ")))
  }
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
	stop("Wrong number of arguments")
}
dir <- args[1]
output_matrix_name <- args[2]

check_files_exist(dir)


z<-ReadMtx(mtx=paste0(dir, "/", "matrix.mtx"), features=paste0(dir, "/", "features.tsv"), cells=paste0(dir, "/", "barcodes.tsv"),feature.column=1)
seurat_object <- CreateSeuratObject(counts = z, project = "Lybestky_project")
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
seurat_object <- ScaleData(seurat_object, features = rownames(seurat_object))
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
pcs_to_use <- 10
# Run UMAP
seurat_object <- RunUMAP(seurat_object, dims = 1:pcs_to_use)
seurat_object <- FindNeighbors(seurat_object, dims = 1:10)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)
all.markers <- FindAllMarkers(seurat_object)

# Connect to Ensembl BioMart
ensembl <- useEnsembl(biomart = "genes")

# Set dataset for pig
ensembl_pig <- useEnsembl(biomart = "genes", dataset = "sscrofa_gene_ensembl")

# Get pig genes from your Seurat object
pig_genes <- rownames(seurat_object)

# Retrieve pig genes with their human orthologs
gene_mapping <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name",
                 "hsapiens_homolog_ensembl_gene", "hsapiens_homolog_associated_gene_name"),
  filters = "external_gene_name",
  values = pig_genes,
  mart = ensembl_pig
)

# Rename columns for clarity
colnames(gene_mapping) <- c("pig_ensembl_gene_id", "pig_gene_name",
                            "human_ensembl_gene_id", "human_gene_name")

# Remove entries without human orthologs
gene_mapping <- gene_mapping[gene_mapping$human_gene_name != "", ]

# Remove duplicate mappings based on pig gene names
gene_mapping <- gene_mapping[!duplicated(gene_mapping$pig_gene_name), ]

# Connect to Ensembl BioMart
ensembl <- useEnsembl(biomart = "genes")

# Set dataset for pig
ensembl_pig <- useEnsembl(biomart = "genes", dataset = "sscrofa_gene_ensembl")

# Get pig genes from your Seurat object
pig_genes <- rownames(seurat_object)

# Retrieve pig genes with their human orthologs
gene_mapping <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name",
                 "hsapiens_homolog_ensembl_gene", "hsapiens_homolog_associated_gene_name"),
  filters = "external_gene_name",
  values = pig_genes,
  mart = ensembl_pig
)

# Rename columns for clarity
colnames(gene_mapping) <- c("pig_ensembl_gene_id", "pig_gene_name",
                            "human_ensembl_gene_id", "human_gene_name")

# Remove entries without human orthologs
gene_mapping <- gene_mapping[gene_mapping$human_gene_name != "", ]

# Remove duplicate mappings based on pig gene names
gene_mapping <- gene_mapping[!duplicated(gene_mapping$pig_gene_name), ]

# Subset Seurat object to genes that have human orthologs
common_genes <- intersect(rownames(seurat_object), gene_mapping$pig_gene_name)
seurat_object <- subset(seurat_object, features = common_genes)

# Create a named vector for mapping pig gene names to human gene names
gene_map <- setNames(gene_mapping$human_gene_name, gene_mapping$pig_gene_name)

# Rename genes in Seurat object
new_gene_names <- gene_map[rownames(seurat_object)]
rownames(seurat_object@assays$RNA@counts) <- new_gene_names
rownames(seurat_object@assays$RNA@data) <- new_gene_names

# Remove any genes with NA after mapping (if any)
valid_genes <- !is.na(rownames(seurat_object))
seurat_object <- subset(seurat_object, features = rownames(seurat_object)[valid_genes])

# Ensure clusters are defined
# If not already clustered, perform clustering
# seurat_object <- FindNeighbors(seurat_object, dims = 1:10)
# seurat_object <- FindClusters(seurat_object, resolution = 0.5)

# Find marker genes for each cluster
markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# View the markers
head(markers)
# Prepare marker gene list for scCATCH
# scCATCH requires a list where each element is a data frame of markers for a cluster
marker_list <- split(markers, markers$cluster)
# Create scCATCH input object
# scCATCH uses the normalized expression data
expr_matrix <- as.matrix(seurat_object@assays$RNA@data)

# Initialize scCATCH object
scCATCH_obj <- createscCATCH(data = expr_matrix, cluster = seurat_object@meta.data$seurat_clusters)

# Add marker genes to scCATCH object
scCATCH_obj@markers <- marker_list
# Set species and tissue for annotation
species <- "Human"
tissue <- "Liver"
# Identify cell types
scCATCH_obj <- findmarkergene(object = scCATCH_obj,
                              species = species,
                              marker = cellmatch,
                              tissue = tissue,
                              use_method = "median",
                              marker_list = NULL)

# View the cell type annotations
celltype_results <- scCATCH_obj@celltype
# Create a data frame mapping clusters to predicted cell types
cluster_annotations <- data.frame(cluster = celltype_results$cluster,
                                  cell_type = celltype_results$cell_type,
                                  stringsAsFactors = FALSE)

# Add cell type annotations to Seurat object metadata
seurat_object@meta.data$scCATCH_cell_type <- cluster_annotations$cell_type[match(seurat_object@meta.data$seurat_clusters, cluster_annotations$cluster)]
