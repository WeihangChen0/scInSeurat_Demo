
suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
  library(Matrix)
  library(tidyverse)
  library(SCENIC)
  library(SCopeLoomR)
  library(Seurat)
  library(gtools)
  library(foreach)
  library(rjson)
})

# cellInfo <- read.csv(cellAnnotation_file)

getpyScenicResultToR <- function(adjacency_matrix_file, 
                                 TF2regulon_file, 
                                 auc_mtx_file){
  ################################################################# 
  ## load genie3/GRNBoost result adjacency_matrix_file
  if (ncol(fread(adjacency_matrix_file))==3){
    adj_mtx <- fread(adjacency_matrix_file)
    names(adj_mtx) <- c("TF","target","importance")
  }
  #################################################################
  
  
  #################################################################                                
  ## load TF2regulon_file
  tmp <- fread(TF2regulon_file)
  
  # split header rows from data
  header_info <- tmp[1:3,]
  TF2enrichment_df <- tmp[-c(1:3),]
  
  # check if column 1 is TF
  if (grepl("TF|transcription_factor",header_info[3,1],ignore.case = T)){
    message('Found transcription factor column')
    colnames(TF2enrichment_df)[1] <- 'TF'
  }
  # check if column 2 is motif
  if (grepl("motif",header_info[3,2],ignore.case = T)){
    message('Found motif_name/motif_ID')
    colnames(TF2enrichment_df)[2] <- 'MotifID'
  }
  
  header1<-setdiff(unique(unlist(header_info[1,])), "")
  header2<-setdiff(unique(unlist(header_info[2,])), "")
  
  # check if row 1 is text "Enrichment"
  if (all(grepl("enrichment",header1,ignore.case = T))){        
    message("line 1 of regulon file is correct")
    
    # check if there is any missing column name on row 2
    if (length(header2)==(ncol(TF2enrichment_df)-2)){
      message("Number of columns of regulon file is correct")
      colnames(TF2enrichment_df)[-c(1:2)] <- paste0("Enrichment_", header2)
    }
  }
  rm(tmp)
  #################################################################  
  
  #################################################################  
  ## load AUCell matrix (and others)
  
  if (grepl('loom$',auc_mtx_file)){
    lfile <- open_loom(auc_mtx_file)
    dgem <- get_dgem(lfile)
    regulons <- get_regulons(lfile, column.attr.name = "Regulons")
    
    regulons_AUC <- get_regulons_AUC(lfile,
                                     column.attr.name = "RegulonsAUC",
                                     rows = "regulons", columns = "Cell"
    )
    
    regulonsAucThresholds <- get_regulon_thresholds(lfile)
    
    close_loom(lfile)
    
  } else if (grepl('csv$',auc_mtx_file)){
    # check if the csv file has been transposed 
    # if so the cell barcodes are used as columns
    tmp <- fread(auc_mtx_file)
    if (sum(grepl("[ATCG]{16}",colnames(tmp)))==(ncol(tmp)-1) ){
      # the output matrix has cell/barcode as column and regulon as row
      regulons_AUC <- as.matrix(t(tmp[,-1]))
      colnames(regulons_AUC) <- as.character(unlist(tmp[,1]))
      rm(tmp)
    } else if (sum(grepl('[ATCG]{16}',unlist(tmp[,1])))==nrow(tmp) ) {
      # the first column is the barcode
      regulons_AUC <- as.matrix(tmp[,-1])
      row.names(regulons_AUC) <- as.character(unlist(tmp[,1]))
      rm(tmp)
    } else {
      regulons_AUC <- NULL
      message(paste0(auc_mtx_file, 
                     " format not consistent with current version"))
    }
    
  }
  
  if (is.null(regulons) | is.null(regulonsAucThresholds)){
    # read csv AUCell 
    outputList <- list(adj_mtx, TF2enrichment_df, regulons_AUC)
    names(outputList)
  } else {
    # read loom AUCell 
    outputList <- list(adj_mtx, TF2enrichment_df, regulons_AUC,
                       regulons, regulonsAucThresholds)
  }
  return(outputList)
}

regulonOnOFF_per_celltype <- function(cellInfo=cellInfo,
                                      cellsAUC=regulons_AUC,
                                      thresholds=regulonsAucThresholds){
  if (typeof(regulons_AUC)=="S4"){
    DR1 <- t(regulons_AUC@assays@data$AUC)
  } else {
    DR1 <- as.matrix(regulons_AUC)
  }
  # print("###### printing regulons_AUC #######")
  # print(DR1)
  # print("####################################")

  ## get threshold AUC into dataframe
  thrd_df <- cbind(names(regulonsAucThresholds), 
                         regulonsAucThresholds) %>% data.frame()

  for (j in 1:ncol(thrd_df)){    
    if (sum(as.numeric(unlist(thrd_df[,j])), na.rm = T)==0){
      thrd_df[,j] <- as.character(thrd_df[,j])
    } else{
      thrd_df[,j] <- as.numeric(thrd_df[,j])
    }
  }

  ## detect whether threshold value was stored as elements or names of element
  if (all(sapply(thrd_df, class)==c("numeric","character"))){
    colnames(thrd_df) <- c("threshold_AUC", "regulon")
  } else if (
    all(sapply(thrd_df, class)==c("character","numeric"))
  ){
    colnames(thrd_df) <- c("regulon","threshold_AUC")
  } else {
    message("check input regulonsAucThresholds")
  }

  if (!is.null(DR1) && ncol(DR1)==nrow(thrd_df)){
    # match column order of DR1 to row order of thrd_df
    thrd_df <- thrd_df[match(colnames(DR1),thrd_df$regulon),]
    
    # check if each regulon is on/off in each cell 
    tmp <- matrix(,nrow=nrow(DR1),ncol(DR1))
    row.names(tmp) <-row.names(DR1)
    colnames(tmp) <- colnames(DR1)    
    for (j in 1:ncol(tmp)){
      tmp[,j] <- as.numeric(DR1[,j] > thrd_df$threshold_AUC[j])
    }
    # print("###### printing tmp #######")
    # print(DR1)
    # print("####################################")

    # group on/off by celltype
    output_df <- matrix(0, nrow = ncol(DR1),ncol = length(unique(cellInfo[,1])))%>%data.frame()
    colnames(output_df) <- mixedsort(as.character(unique(cellInfo[,1])))
    row.names(output_df)<- colnames(DR1)
    
    for (j in 1:ncol(output_df)){
      b <- row.names(cellInfo)[cellInfo[,1]==colnames(output_df)[j]]
      if (length(b)>0){
        pct_on <- colSums(tmp[row.names(tmp)%in%b,])/length(b)
        pct_on[is.na(pct_on)] <- 0
        output_df[,j] <- round(pct_on,3)

      } else {
        output_df[,j] <- 0
      }
    }
  }
  return(output_df)
   
}

#* This function is to return a table of gene and its corresponding weight/importance
#* in a given regulon. The user need to specify the regulon name and supply the
#* output of pyscenic ctx command
getGeneWeightInRegulon <- function(TF2enrichment_df,
                                   regulon_name){
  
  #* if regulon name contains special string to 
  #* indicate extended regulon, strip it off from name 
  regex_pat <- "\\[|\\(|\\)|\\]|\\+"
  if(grepl(regex_pat, regulon_name)){
    regulon_name <- gsub(regex_pat,"",regulon_name)
  }
  
  # get target gene and corresponding importance in list
  TargetGenes <- TF2enrichment_df$Enrichment_TargetGenes[TF2enrichment_df$TF==regulon_name]
  regulon_enrich_list<- strsplit(TargetGenes, split="), (", fixed = T)
  
  if (length(regulon_enrich_list)==0){
    message("regulon not found")
  } else {
    
    regulon_enrich_df <- foreach(i=1:length(regulon_enrich_list), 
                                 .combine=rbind)%do%{

      tmp <- data.frame(regulon_enrich_list[[i]]) %>% 
        dplyr::rename(X1 = 1)%>% 
        separate("X1",c('gene','importance'),sep="', ")%>% 
        mutate(gene_name = gsub("'", "" , gene, fixed = T) )%>%
        select(gene_name, importance)
      
      tmp$set<-i
      
      return(tmp)
    }
    regulon_enrich_df$gene_name <- gsub(regex_pat,"",
                                        regulon_enrich_df$gene_name)
    regulon_enrich_df$importance <- gsub(regex_pat,"",
                                         regulon_enrich_df$importance)
    regulon_enrich_df$importance <- signif(as.numeric(regulon_enrich_df$importance),3)
    
    output_by_regulon <- unique(regulon_enrich_df[,c('gene_name','importance')])%>%
      arrange(desc(importance))
  }
  
  return(output_by_regulon)
  
}




geneRelatedRegulonInClusters<-function(cellInfo,TF2enrichment_df,
                  regulons_AUC,regulonsAucThresholds,
                  geneName){
  
  regulon_score_celltype<-regulonOnOFF_per_celltype(cellInfo=cellInfo,
                                                    cellsAUC=regulons_AUC,
                                                    thresholds=regulonsAucThresholds)
  
  # get regulons that have the gene
  rgs <- unique(TF2enrichment_df$TF[grep(geneName,
                                         TF2enrichment_df$Enrichment_TargetGenes)])
  
  geneWeightInRegulons<- foreach(rg=rgs, .combine = rbind)%do%{
    tmp <- getGeneWeightInRegulon(TF2enrichment_df=TF2enrichment_df,
                                  regulon_name = rg)
    
    tmp$regulon_name <- rg
    tmp$rank <- which(tmp$gene_name==geneName)
    tmp$regulon_size <- (nrow(tmp)-1)
    
    return(tmp[tmp$gene_name==geneName,])
  }
  geneWeightInRegulons <- geneWeightInRegulons %>% 
    mutate(i2 = rank/regulon_size)%>%
    arrange(i2)%>%select(-i2)
  
  geneWeightInRegulons$highlyAffectedCluster <- ""

  for (i in 1:nrow(geneWeightInRegulons)){
    rgn <- geneWeightInRegulons$regulon_name[i]
    
    # get activity of all cluster of this regulon
    ind <- which(startsWith(row.names(regulon_score_celltype), rgn))
    
    if (length(ind)==1){
      tmp <- sort(unlist(regulon_score_celltype[ind,]), decreasing = T)
      ## output top 5 clusters only
      # if (length(tmp>5)){
      #   tmp <- tmp[1:5]
      # }
      geneWeightInRegulons$highlyAffectedCluster[i] <- rjson::toJSON(tmp)
    } else if (length(ind)!=1){
      message("regulon not found")
    }

  }
  return(geneWeightInRegulons)
}
