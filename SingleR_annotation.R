# Install packages if not already installed
if (!require("Seurat")) install.packages("Seurat")
if (!require("SingleR")) BiocManager::install("SingleR")
if (!require("celldex")) BiocManager::install("celldex")
if (!require("biomaRt")) install.packages("biomaRt")
if (!require("dplyr")) install.packages("dplyr")

# Load libraries
library(Seurat)
library(SingleR)
library(celldex)
library(biomaRt)
library(dplyr)

library(Seurat)
z<-ReadMtx(mtx="matrix.mtx", features="features.tsv", cells="barcodes.tsv",feature.column=1)
seurat_object <- CreateSeuratObject(counts = z, project = "HumanLiver")
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

# Load your Seurat object (replace with your actual file)
# load("your_seurat_object.RData")
# Alternatively, if saved as an RDS file
# your_seurat_object <- readRDS("your_seurat_object.rds")

# Use biomaRt to map genes
# Connect to Ensembl BioMart
pig_mart <- useEnsembl("ensembl", dataset = "sscrofa_gene_ensembl")
human_mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")

# Get mapping between pig and human genes
gene_mapping <- getLDS(attributes = c("ensembl_gene_id", "external_gene_name"),
                       filters = "ensembl_gene_id",
                       values = rownames(seurat_object),
                       mart = pig_mart,
                       attributesL = c("ensembl_gene_id", "external_gene_name"),
                       martL = human_mart,
                       uniqueRows = TRUE)

# Rename columns for clarity
colnames(gene_mapping) <- c("pig_ensembl_gene_id", "pig_gene_name",
                            "human_ensembl_gene_id", "human_gene_name")

# Remove duplicate mappings
gene_mapping <- gene_mapping[!duplicated(gene_mapping$pig_gene_name), ]

# Subset Seurat object to genes that have human orthologs
common_genes <- intersect(rownames(seurat_object), gene_mapping$pig_gene_name)
seurat_object <- subset(seurat_object, features = common_genes)

# Update gene names to human orthologs
# Create a named vector for mapping
gene_map <- setNames(gene_mapping$human_gene_name, gene_mapping$pig_gene_name)
# Rename genes in Seurat object
rownames(seurat_object@assays$RNA@counts) <- gene_map[rownames(seurat_object@assays$RNA@counts)]
rownames(seurat_object@assays$RNA@data) <- gene_map[rownames(seurat_object@assays$RNA@data)]

# Load Human Primary Cell Atlas Data
hpca.se <- celldex::HumanPrimaryCellAtlasData()
# Extract the normalized data matrix
data_matrix <- GetAssayData(seurat_object, slot = "data")

# Run SingleR
singleR_results <- SingleR(test = data_matrix,
                           ref = hpca.se,
                           labels = hpca.se$label.main,
                           de.method = "wilcox")

# View results
head(singleR_results$labels)
# Add SingleR labels to Seurat metadata
seurat_object <- AddMetaData(seurat_object, metadata = singleR_results$labels, col.name = "SingleR.labels")

# If you have clusters and want to annotate clusters
cluster_labels <- data.frame(cluster = levels(Idents(seurat_object)))
cluster_labels$SingleR.labels <- singleR_results$pruned.labels
# Map cluster annotations back to Seurat object
seurat_object@meta.data$SingleR.cluster.labels <- cluster_labels$SingleR.labels[match(Idents(seurat_object), cluster_labels$cluster)]
# Save the updated Seurat object
saveRDS(seurat_object, file = "seurat_object_annotated.rds")
