library(tidyverse)
library(Seurat)

# load seurat object
allSample_BC <- readRDS('tmp_result/Seurat_not_integrated.rds')

# for res=0.3, do marker calling for all clusters 
Idents(allSample_BC) <- 'RNA_snn_res.0.3'
combined_marker<-FindAllMarkers(allSample_BC)
combined_marker$avg_log2FC <- round(combined_marker$avg_log2FC,3)
combined_marker$p_val <- signif(combined_marker$p_val, 4)
combined_marker$p_val_adj <- signif(combined_marker$p_val_adj, 4)

write.csv(combined_marker,
          'Result/allMarker.res0.3.csv',
          quote = F,row.names = F)

