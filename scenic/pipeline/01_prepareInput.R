suppressPackageStartupMessages({
  library(readxl)
  library(Matrix)
  library(tidyverse)
  library(SCENIC)
  library(SCopeLoomR)
  library(Seurat)
})
test_data_dir <- "~/test_scenic_metabolic" 
# load raw
seuratObj <- readRDS(file.path(test_data_dir, "raw_data/2021-02-03_seuratObj.Rds"))
cts <- seuratObj@assays$RNA@counts
metadata <- seuratObj@meta.data
rm(seuratObj)


# load target genes
metabolic_genes <- read_excel(file.path(test_data_dir,"TFs_and_metabolic_enzymes_for scenic_analysis.xlsx"))
allTFs_dmel <- read_csv(file.path(test_data_dir,"resources/TFs/allTFs_dmel.txt"),
                        col_names = F)

# check gene names
setdiff(metabolic_genes$Gene, row.names(cts))
# gene "Lime" was "CG18446" in the previous version, use this name instead 
metabolic_only_genes <- c(intersect(metabolic_genes$Gene, row.names(cts)),
                                    "CG18446")
write_tsv(data.frame(intersect(allTFs_dmel$X1, metabolic_only_genes)),
          "resources/TFs/metabolic_TFs_dmel.txt", col_names=F)

# build loom file with metabolic genes and selected highly variable genes
build_loom(file.path(test_data_dir,"Data/8W_exprMat_regular.loom"),
           cts[gene_keep,])
close_loom(file.path(test_data_dir,"Data/8W_exprMat_regular.loom"))

