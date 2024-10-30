# Install packages if not already installed
if (!require("Seurat")) install.packages("Seurat")
if (!require("SingleR")) BiocManager::install("SingleR")
if (!require("celldex")) BiocManager::install("celldex")
if (!require("biomaRt")) install.packages("biomaRt")
if (!require("dplyr")) install.packages("dplyr")
if (!require("anndata")) BiocManager::install("anndata")


library(Seurat)
library(SingleR)
library(celldex)
library(biomaRt)
library(dplyr)
library(anndata)

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
if (length(args) != 3) {
	stop("Wrong number of arguments")
}
dir <- args[1]
output_matrix_name <- args[3]

check_files_exist(dir)

specie <- args[2]


z<-ReadMtx(mtx=paste0(dir, "/", "matrix.mtx"), features=paste0(dir, "/", "features.tsv"), cells=paste0(dir, "/", "barcodes.tsv"),feature.column=1)

specie <- "sscrofa"
ensembl <- useEnsembl(biomart = "genes")

# Set datasets for pig and human

ensembl_pig <- useEnsembl(biomart = "genes", dataset = paste0(specie, "_gene_ensembl"), host="https://may2024.archive.ensembl.org")
ensembl_human <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", host="https://may2024.archive.ensembl.org")

pig_genes <- rownames(z)
gene_mapping <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name",
                 "hsapiens_homolog_ensembl_gene", "hsapiens_homolog_associated_gene_name"),
  filters = "ensembl_gene_id",
  values = pig_genes,
  mart = ensembl_pig
)

colnames(gene_mapping) <- c("specie_ensembl_gene_id", "specie_gene_name",
                            "human_ensembl_gene_id", "human_gene_name")
gene_mapping <- gene_mapping[gene_mapping$human_gene_name != "", ]
gene_mapping <- gene_mapping[!duplicated(gene_mapping$specie_gene_name), ]

gene_map <- setNames(gene_mapping$human_gene_name, gene_mapping$specie_ensembl_gene_id)
new_gene_names <- gene_map[rownames(z)]
valid_genes <- !is.na(new_gene_names)
z <- z[valid_genes, ]

seurat_object <- CreateSeuratObject(counts = z, project = "HumanLiver")
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
seurat_object <- ScaleData(seurat_object, features = rownames(seurat_object))
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
hpca.se <- celldex::HumanPrimaryCellAtlasData()
new_gene_names <- new_gene_names[valid_genes]
rownames(z) <- new_gene_names
singleR_results <- SingleR(test = z,
                           ref = hpca.se,
                           labels = hpca.se$label.main,
                           de.method = "wilcox")
labels <- as.factor(singleR_results$labels)
names(labels) <- rownames(seurat_object@assays$RNA@cells)
seurat_object@active.ident <- labels
all.markers <- FindAllMarkers(seurat_object)
write.table(x=all.markers, file=output_matrix_name)
