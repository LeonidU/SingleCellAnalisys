library(Seurat)
library(anndata)


args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
z<-ReadMtx(mtx=paste0(dir, "/", "matrix.mtx"), features=paste0(dir, "/", "features.tsv"), cells=paste0(dir, "/", "barcodes.tsv"),feature.column=1)
seurat_object <- CreateSeuratObject(counts = z, project = "Lybestky_project")
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
mat <- as.matrix(seurat_object@assays$RNA$data)
seurat_object <- ScaleData(seurat_object, features = rownames(seurat_object))
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
pcs_to_use <- 10

seurat_object <- RunUMAP(seurat_object, dims = 1:pcs_to_use)
seurat_object <- FindNeighbors(seurat_object, dims = 1:10)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)
all.markers <- FindAllMarkers(seurat_object)

x <- AnnData(X=mat, var=as.data.frame(seurat_object@active.ident))
output_matrix_name <- paste0(dir,".h5ad")
write_h5ad(x, output_matrix_name)
write.table(x = all.markers, file = args[2])
