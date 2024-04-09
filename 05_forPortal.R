library(tidyverse)
library(Seurat)
library(DropletUtils)


source('~/Desktop/forSeurat/reformatBarcode.r')


seuratObj <- readRDS('tmp_result/Seurat_not_integrated.rds')

# export umap and metadata 
umap <- seuratObj@reductions$umap@cell.embeddings
metadata <- seuratObj@meta.data

# rename barcodes for exported files 
barcodes <- reformatBarcode(colnames(seuratObj))

row.names(metadata) <- barcodes
row.names(umap) <- barcodes
write.csv(metadata, "Result/metadata.csv")
write.csv(umap, "Result/umap.csv")




