library(DropletUtils)
library(Seurat)
library(Matrix)
library(Seurat)
library(cowplot)
library(SeuratDisk)
library(tidyverse)
library(zellkonverter)
library(reticulate)

## load genome annotation 
dmel_all_r6_56 <- read_delim("dmel-all-r6.56.gtf",col_names = F,delim="\t")%>%
  filter(X3=='gene')

## get gene name to gene id in dataframe
FBgn2Name <-dmel_all_r6_56%>%filter(X3=='gene')%>%
  separate(X9, c('gene_id','gene_name'),sep='; ')%>%
  select(gene_id,gene_name)
FBgn2Name$gene_id <- gsub('gene_id ', '', gsub('"', '', FBgn2Name$gene_id, fixed = T))
FBgn2Name$gene_name <- gsub('gene_symbol ', '', gsub('"', '', FBgn2Name$gene_name, fixed = T))
FBgn2Name$gene_name <-gsub(';', '', FBgn2Name$gene_name, fixed = T)
rm(dmel_all_r6_56)


########################################################################
# load all samples in Data/ into seurat object
samples <- list.dirs("Data",recursive = F,full.names = F)
seurat_list <- list()

for (i in 1:length(samples)){
  # load all barcode passed whitelist
  barcodes <- read_csv(file.path("Data",
                                 samples[i],
                                 "Solo.out/GeneFull/raw/barcodes.tsv"),
                       col_names = F)
  # load barcode filtered for empty 
  barcodes_filtered <- read_csv(file.path("Data",
                                          samples[i],
                                          "Solo.out/GeneFull/filtered/barcodes.tsv"),
                       col_names = F)
  
  # load features
  features <- read_delim(file.path("Data",
                                   samples[i],
                                   "Solo.out/GeneFull/raw/features.tsv"), 
                         delim = "\t",col_names = F)
  features$X2 <- FBgn2Name$gene_name[match(features$X1,FBgn2Name$gene_id)]
  
  # load h5ad  
  ad <- import("anndata", convert = FALSE)
  test_ad <- ad$read_h5ad(file.path("Data",
                                    samples[i],
                                    "Solo.out/GeneFull/raw/UniqueAndMult.h5ad"))
  sce_test <- AnnData2SCE(test_ad)
  # remove index 0 row/column
  x <- t(sce_test@assays@data$X[-1,-1])
  row.names(x) <-features$X2
  colnames(x) <- barcodes$X1

  # save to seurat  
  seurat_list[[samples[i]]] <- CreateSeuratObject(x[,which(colnames(x) %in% barcodes_filtered$X1)])
}

tmp <- merge(seurat_list[[1]],seurat_list[2:length(samples)])

# remove low expressed gene
allSample_raw<- CreateSeuratObject(tmp@assays$RNA@counts[
  rowSums(tmp@assays$RNA@counts)>10,])
rm(tmp)

saveRDS(allSample_raw, "tmp_result/allSample_raw.rds")
