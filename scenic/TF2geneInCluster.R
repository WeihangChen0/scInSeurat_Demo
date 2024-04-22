cellInfo <- read.csv("cellInfo.W.csv", row.names=1)

outputList <- getpyScenicResultToR(adjacency_matrix_file="Output/W.regular-metabolite.adj.tsv",
                     TF2regulon_file="Output/W.regulons.tsv",
                     auc_mtx_file="Output/8W.auc_mtx.loom")
adj_mtx <- outputList[[1]]
TF2enrichment_df <- outputList[[2]]
regulons_AUC <- outputList[[3]]
if (length(outputList)==5){
  regulons <- outputList[[4]]
  regulonsAucThresholds <- outputList[[5]]
}
rm(outputList)
# 
tmp <- getGeneWeightInRegulon(TF2enrichment_df,"B-H1(+)")
regulon_score_celltype<-regulonOnOFF_per_celltype(cellInfo=cellInfo,
                                                  cellsAUC=regulons_AUC,
                                                  thresholds=regulonsAucThresholds)
