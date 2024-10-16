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

# Retrieve pig genes with their human orthologs
# Note: We can retrieve homologs directly using the pig dataset

ensembl <- useEnsembl(biomart = "genes")

# List available datasets (optional)
# listDatasets(ensembl)

# Set datasets for pig and human
ensembl_pig <- useEnsembl(biomart = "genes", dataset = "sscrofa_gene_ensembl")
ensembl_human <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Get pig genes from your Seurat object
pig_genes <- rownames(seurat_object)

# Retrieve pig genes with their human orthologs
# Note: We can retrieve homologs directly using the pig dataset
gene_mapping <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name",
                 "hsapiens_homolog_ensembl_gene", "hsapiens_homolog_associated_gene_name"),
  filters = "ensembl_gene_id",
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
#gene_map <- setNames(gene_mapping$human_gene_name, gene_mapping$pig_gene_name)
 gene_map <- setNames(gene_mapping$human_gene_name, gene_mapping$pig_ensembl_gene_id)

counts_matrix <- GetAssayData(seurat_object, assay = "RNA", slot = "counts")

# Map the row names (gene names) to human orthologs
new_gene_names <- gene_map[rownames(counts_matrix)]

# Remove any genes where the mapping resulted in NA
valid_genes <- !is.na(new_gene_names)
counts_matrix <- counts_matrix[valid_genes, ]
new_gene_names <- new_gene_names[valid_genes]

# Update row names with human gene names
#rownames(counts_matrix) <- new_gene_names

# Set the modified counts matrix back into the Seurat object
#seurat_object <- SetAssayData(
#  object = seurat_object,
#  slot = "counts",
#  new.data = counts_matrix
#)

# Repeat the process for the "data" slot (normalized data)
data_matrix <- GetAssayData(seurat_object, assay = "RNA", slot = "scale.data")

# Subset and rename the data matrix
data_matrix <- data_matrix[valid_genes, ]
rownames(data_matrix) <- new_gene_names

# Set the modified data matrix back into the Seurat object
#seurat_object <- SetAssayData(
#  object = seurat_object,
#  slot = "data",
#  new.data = data_matrix
#)

# Load Human Primary Cell Atlas Data
hpca.se <- celldex::HumanPrimaryCellAtlasData()
# Extract the normalized data matrix
#data_matrix <- GetAssayData(seurat_object, slot = "data")

# Run SingleR
singleR_results <- SingleR(test = data_matrix,
                           ref = hpca.se,
                           labels = hpca.se$label.main,
                           de.method = "wilcox")

# View results
#head(singleR_results$labels)
# Add SingleR labels to Seurat metadata
seurat_object <- AddMetaData(seurat_object, metadata = singleR_results$labels, col.name = "SingleR.labels")

# If you have clusters and want to annotate clusters
#cluster_labels <- data.frame(cluster = levels(Idents(seurat_object)))
#cluster_labels$SingleR.labels <- singleR_results$pruned.labels
# Map cluster annotations back to Seurat object
#seurat_object@meta.data$SingleR.cluster.labels <- cluster_labels$SingleR.labels[match(Idents(seurat_object), cluster_labels$cluster)]
# Save the updated Seurat object
#saveRDS(seurat_object, file = "seurat_object_annotated.rds")

library(anndata)
new_list <- list()
#It's kind of magic, but without that loading into h5ad doesn't work
for (i in seq(1, length(singleR_results$pruned.labels))) {
	new_list[as.character(i)] <- singleR_results$pruned.labels[i]
}
mat <- as.matrix(data_matrix)
rownames(mat)[duplicated(rownames(mat))] <- paste0(rownames(mat)[duplicated(rownames(mat))], ".1")

x <- AnnData(X=mat, uns=new_list, var=as.data.frame(seurat_object@active.ident))

write_h5ad(x, "pigs.h5ad")
