library(tidyverse)
library(gridExtra)
library(Cairo)
library(Seurat)

# set input file path 
adjacency_matrix_file <- ""
TF2regulon_file <- ""
auc_mtx_file <- ""
cellAnnotation_file <- ""


# cellInfo <- read.csv(cellAnnotation_file)

# load scenic result 
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


# calculate and plot RSS for all regulons per celltype
rssPlotPerCelltype <- function(regulons_AUC, cellInfo,
                               output_format="svg"){
    # calculate rss per celltype
    rss <- calcRSS(AUC = regulons_AUC,
                   cellAnnotation = cellInfo, cellTypes =cellInfo)
    rss <- rss[,gtools::mixedsort(colnames(rss))]

    # plot 1 rss plot for each cell type  
    plot_list <- list()
    for (i in 1:ncol(rss)){
        # make data frame for ggplot
        pdata <- cbind(sort(rss[,i],decreasing=T),
                       1:nrow(rss)) %>% data.frame()
        pdata <- cbind(row.names(pdata),pdata)
        names(pdata) <- c('regulon', 'rss',"order")
        pdata$cluster <- colnames(rss)[i]
        pdata <- pdata %>% filter(rss>0.01) # remove less important regulons to this cluster
        
        # dotplot of curve for rss and label the top 6 regulon with text
        p <- ggplot(pdata, aes(x= order, y= rss, 
                               colour="green",label=regulon))+
                geom_point() +
                geom_label_repel(aes(label=ifelse(rss>pdata[6,'rss'],
                                                    as.character(regulon),'')),
                                hjust=0,vjust=0) +
                theme(axis.title.x=element_blank(),
                        axis.text.x=element_blank(),
                        axis.ticks.x=element_blank())+
                labs(title = colnames(rss)[1], y = "RSS")
        
        plot_list[[i]] <- p
    }
    # write plots to file
    if (output_format=="svg"){
        svg(file.path('tmp_result',"rss_per_celltype.svg"))
        do.call("grid.arrange", c(plot_list, ncol = 1))
        dev.off()
    } else if (output_format=="pdf"){
        pdf(file.path('tmp_result',"rss_per_celltype.pdf"))
        do.call("grid.arrange", c(plot_list, ncol = 1))
        dev.off()
    } else if (output_format=="png"){
        png(file.path('tmp_result',"rss_per_celltype.png"))
        do.call("grid.arrange", c(plot_list, ncol = 1))
        dev.off()
    } else {
       message("this function export plot in pdf/svg/png format only")
    }
}

# do UMAP based on regulon and plot regulon activity on this umap
AUC_UMAP <- function(regulons_AUC, exprMat = dgem, 
                     cellsAUC=regulons_AUC,
                     thresholds=regulonsAucThresholds,
                     plot_type = c("binaryAUC", "AUC"),
                     output_pdf_file){

    if (typeof(regulons_AUC)=="S4"){
        try(DR1 <- t(regulons_AUC@assays@data$AUC))
    } else {
        DR1 <- as.matrix(regulons_AUC)
    }

    scenic_umap <- RunUMAP(DR1,
                           umap.method = "umap-learn",
                           n.neighbors=10, min_dist=0.4)

    CairoPDF(output_pdf_file, width=20, height=15)
    par(mfrow=c(4,6))
    AUCell_plotTSNE(tSNE=scenic_umap@cell.embeddings,
                    exprMat = dgem, plots=plot_type,
                    cellsAUC=DR1,
                    thresholds=regulonsAucThresholds)
    dev.off()
}
