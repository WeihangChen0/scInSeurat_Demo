source('~/Desktop/forSeurat/DEGcallingPerCluster.R')

# load seurat object with known clustering
allSample_BC <- readRDS('tmp_result/Seurat_not_integrated.rds')

# compare treated vs untreated/widltype group in all clusters
output <- DEGcallingPerCluster(allSample_BC,
                               control_group=c("wt"),
                               group_name="condition",
                               idents_name= 'RNA_snn_res.0.3',n_cores=3)

write_tsv(output,
          'Result/DEG_all_condition.res0.3.tsv')