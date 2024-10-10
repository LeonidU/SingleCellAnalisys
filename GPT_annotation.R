# GPT4 cell types annotation
#
# Command line format
# Rscript GPT_annotation.R raw_data_folder OpenAI_API_KEY output_matrix_name
# raw_data_folder - folder with cellranger output, must contain matrix.mtx, features.tsv, barcodes.tsv
# output_matrix_name - h5ad filename for output

library(Seurat)
library(anndata)
# Load packages
library(GPTCelltype)
library(openai)

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
openai_key <- args[2]
output_matrix_name <- args[3]

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


Sys.setenv(OPENAI_API_KEY = openai_key)
# Assume you have already run the Seurat pipeline https://satijalab.org/seurat/
# "obj" is the Seurat object; "markers" is the output from FindAllMarkers(obj)
# Cell type annotation by GPT-4


res <- gptcelltype(all.markers, model = 'gpt-4o')

# Assign cell type annotation back to Seurat object
seurat_object@meta.data$celltype <- as.factor(res[as.character(Idents(seurat_object))])
mat <- as.matrix(seurat_object@assays$RNA$scale.data)
x <- AnnData(X=mat, uns=as.list(as.data.frame(seurat_object@meta.data$celltype)), var=as.data.frame(seurat_object@active.ident))
write_h5ad(x, output_matrix_name)

