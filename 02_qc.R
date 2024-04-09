# 01 is run soupX on server.

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(harmony)
  library(SCopeLoomR)
  library(patchwork)
  library(SoupX)
  library(Signac)
  library(gridExtra)
})

in_data_dir <- './Data'
samples <- list.dirs(in_data_dir, recursive=F, full.names = F)
sam_grp <- c('wt','wt','foxo','foxo','reptor','reptor')


# in_data_dir <- 'Data'
seurat_list <- lapply(1:6, function(i){
  smpl<-samples[i]
  cur_data <- Read10X(file.path(in_data_dir,smpl,'outs/filtered_feature_bc_matrix'))
  cur_seurat <- CreateSeuratObject(
    counts = cur_data,
    min.cells=3
  )
  cur_seurat$SampleID <- smpl
  
  # Novelty score
  cur_seurat$log10GenesPerUMI <- log10(cur_seurat$nFeature_RNA) / log10(cur_seurat$nCount_RNA)
  cur_seurat$condition <- sam_grp[i]
  print("---- Plot per dataset before QC and batch correction ----")
  # library size
  p1<- data.frame(cur_seurat@meta.data) %>% 
    ggplot(aes(x=condition, fill=condition)) + 
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle("NCells")
  
  # plot UMI counts (transcripts) per cell
  p2 <- data.frame(cur_seurat@meta.data) %>% 
    ggplot(aes(color=condition, fill= condition, x=nCount_RNA)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density") +
    geom_vline(xintercept = 500)
  
  # plot Genes detected per cell
  p3 <- data.frame(cur_seurat@meta.data) %>% 
    ggplot(aes(color=condition, fill= condition, x=nFeature_RNA)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density") +
    geom_vline(xintercept = 500)
  
  # Complexity
  p4<- data.frame(cur_seurat@meta.data) %>% 
    ggplot(aes(x=log10GenesPerUMI,color=condition, fill= condition)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    geom_vline(xintercept = 0.8)
    
  p5<- data.frame(cur_seurat@meta.data) %>% 
    ggplot(aes(x=nCount_RNA, y=nFeature_RNA,color=log10GenesPerUMI)) + 
    geom_point() + 
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method=lm) +
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() +
    geom_vline(xintercept = 500) +
    geom_hline(yintercept = 250) +
    facet_wrap(~condition)
  plot_list <- list(p1, p2, p3, p4, p5) 
  svg(here('tmp_result',
           paste0(smpl,".prefilterQC.svg")))
  do.call("grid.arrange", c(plot_list, ncol = 1))  
  dev.off()
  return(cur_seurat)
})


allSample_raw <-  merge(x=seurat_list[[1]], 
                        y = seurat_list[2:length(seurat_list)], 
                        add.cell.ids = sam_grp,  
                        project = "FC_07486")

save(allSample_raw, file="tmp_result/unfiltered_seurat.RData")

# load soupX data
seurat_list <- list()
seurat_list <- lapply(1:6, function(i){
  smpl<-samples[i]
  cur_data <- Read10X(file.path(in_data_dir,
                                smpl,'outs/soupX_pbmc10k_filt/'))
  cur_seurat <- CreateSeuratObject(
    counts = cur_data,
    min.cells=3
  )
  cur_seurat$SampleID <- smpl
  
  # Novelty score
  cur_seurat$log10GenesPerUMI <- log10(cur_seurat$nFeature_RNA) / log10(cur_seurat$nCount_RNA)
  cur_seurat$condition <- sam_grp[i]
  print("---- Plot per dataset before QC and batch correction ----")
  # library size

  p1<- data.frame(cur_seurat@meta.data) %>% 
    ggplot(aes(x=condition, fill=condition)) + 
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle("NCells")
  
  # plot UMI counts (transcripts) per cell
  p2 <- data.frame(cur_seurat@meta.data) %>% 
    ggplot(aes(color=condition, fill= condition, x=nCount_RNA)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density") +
    geom_vline(xintercept = 500)
  
  # plot Genes detected per cell
  p3 <- data.frame(cur_seurat@meta.data) %>% 
    ggplot(aes(color=condition, fill= condition, x=nFeature_RNA)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density") +
    geom_vline(xintercept = 500)
  
  # Complexity
  p4<- data.frame(cur_seurat@meta.data) %>% 
    ggplot(aes(x=log10GenesPerUMI,color=condition, fill= condition)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    geom_vline(xintercept = 0.8)
  
  
  p5<- data.frame(cur_seurat@meta.data) %>% 
    ggplot(aes(x=nCount_RNA, y=nFeature_RNA,color=log10GenesPerUMI)) + 
    geom_point() + 
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method=lm) +
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() +
    geom_vline(xintercept = 500) +
    geom_hline(yintercept = 250) +
    facet_wrap(~condition)
  plot_list <- list(p1, p2, p3, p4, p5) 
  svg(paste0(smpl,".soupX.prefilterQC.svg"))
  do.call("grid.arrange", c(plot_list, ncol = 1))  
  dev.off()
  return(cur_seurat)
})

allSample_soupX <-  merge(x=seurat_list[[1]], y = seurat_list[2:length(seurat_list)], 
                        add.cell.ids = sam_grp,  
                        project = "FC_07486")
save(allSample_soupX, file="tmp_result/unfiltered_soupX.RData")

#* Muscle cells were expected to have less genes
#* per cell. The qc cutoff were lowered to 100 instead of 200.
allSample_filtered <- subset(x = allSample_raw, 
                             subset= (nCount_RNA >= 100) & 
                               (nFeature_RNA >= 100) & 
                               (log10GenesPerUMI > 0.8))

# remove genes expressed in less than 10 cells 
counts <- GetAssayData(object = allSample_filtered, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]
allSample_filtered <- CreateSeuratObject(filtered_counts,
                                         meta.data = allSample_filtered@meta.data)%>% NormalizeData()

# Dimension reduction with lsi. 
allSample_normalized <- RunTFIDF(allSample_normalized)
allSample_normalized <- FindTopFeatures(allSample_normalized, min.cutoff = 20)
allSample_normalized <- RunSVD(allSample_normalized)
allSample_normalized <- RunUMAP(allSample_normalized, dims = 1:50, reduction = 'lsi')

# plot UMAP before BC
UMAP_beforeBC<-DimPlot(allSample_normalized, 
                       split.by = 'SampleID', pt.size = 0.1)
svg("UMAP before batch correction.svg",width = 14)
DimPlot(allSample_normalized, split.by = 'SampleID', pt.size = 0.1)
dev.off()

# run harmony 
allSample_BC <- RunHarmony(object = allSample_normalized, 
                           group.by.vars = "SampleID", 
                           reduction = "lsi", assay.use = "RNA", 
                           project.dim = FALSE)

allSample_BC <- RunUMAP(allSample_BC,
                        dims = 1:30, reduction = "harmony"
)
allSample_BC <- FindNeighbors(allSample_BC,
                              reduction = "harmony",
                              dims = 1:30) %>% 
  FindClusters(res = seq(0.1,1,0.1) )


# choose res=0.3 to combine 2 close muscle clusters under res=0.4
Idents(allSample_BC) <- 'RNA_snn_res.0.3'
UMAP_afterBC <- DimPlot(allSample_BC, split.by ='SampleID', 
                      pt.size = 0.1,label = T)
svg("UMAP after batch correction res=0.3.svg",
    width = 14)
print(UMAP_afterBC)
dev.off()

# save seurat object to tmp_result
saveRDS(allSample_BC, file = 'tmp_result/Seurat_not_integrated.rds')


