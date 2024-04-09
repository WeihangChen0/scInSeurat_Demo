#* main script from 
#* https://cellgeni.github.io/notebooks/html/new-10kPBMC-SoupX.html
#* plot and detail explanation from 
#* https://github.com/constantAmateur/SoupX/blob/master/vignettes/pbmcTutorial.Rmd

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Seurat)
  library(tidyverse)
  library(SoupX)
  library(DropletUtils)
})

# move cellranger count result folders to Data/
in_data_dir <- "Data"
# get sampleID per folder
samples <- list.dirs(in_data_dir, recursive=F, full.names = F)

for (samp in samples){
  # read in raw & filtered h5 into dgCMatrix classes
  filt.matrix <- Read10X_h5(file.path(in_data_dir,samp,"outs/filtered_feature_bc_matrix.h5"),use.names = T)
  raw.matrix  <- Read10X_h5(file.path(in_data_dir,samp,"outs/raw_feature_bc_matrix.h5"),use.names = T)
  
  # make a soup channel 
  soup.channel  <- SoupChannel(raw.matrix, filt.matrix)
  
  #* SoupX requires clusters in order to define marker genes. 
  #* Quickly cluster using Seurat.
  srat  <- CreateSeuratObject(counts = filt.matrix)
  srat    <- SCTransform(srat, verbose = F)
  srat    <- RunPCA(srat, verbose = F)
  srat    <- RunUMAP(srat, dims = 1:30, verbose = F)
  srat    <- FindNeighbors(srat, dims = 1:30, verbose = F)
  srat    <- FindClusters(srat, verbose = T)
  
  # add clusters to soup channel 
  meta    <- srat@meta.data
  umap    <- srat@reductions$umap@cell.embeddings
  soup.channel  <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
  soup.channel  <- setDR(soup.channel, umap)
  
  
  # run soupX
  soup.channel  <- autoEstCont(soup.channel)
 
  head(soup.channel$soupProfile[order(soup.channel$soupProfile$est, 
                                      decreasing = T), ], n = 20)
  
  
  # round figure to integer and write the directory with corrected read counts.
  adj.matrix  <- adjustCounts(soup.channel, roundToInt = T)
  write10xCounts(file.path(in_data_dir,
                           samp, "outs/soupX_pbmc10k_filt"), 
                 adj.matrix)
  save.image(file.path(in_data_dir,samp,
                       "outs/soupX_pbmc10k_filt/soupX.RData"))
}