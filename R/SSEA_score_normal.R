#' Spatial genesets enrichment score for normal spots
#'
#'Predicting cell types in normal regions from spatial transcriptomic data.
#'
#' @param HallMarker Gene markers for normal spots.
#' @param test_gene Top_gene lists.
#' @param population_size Number of genes in the species studied.
#' @param method Statistical methods for p-value correction. c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#'
#' @return A cell type prediction list, containing spot and cluster cell type prediction.
#' @export
#'
#' @examples
#' \dontrun{pred <- SSEA_score_normal(HallMarker,test_gene,20000)}
SSEA_score_normal <- function(HallMarker,test_gene,population_size,method="BH"){
  colnames(HallMarker) <- c("ref_type1","ref_type2","makrer")
  SSEA_list <- list()

  ####CT1 anno
  result_df <- c()
  HallMarker1 <- unique(HallMarker[,c("ref_type1","makrer")])
  gene_stat <- data.frame(table(HallMarker1$ref_type1))
  HallMarker1 <- subset(HallMarker1,ref_type1 %in% gene_stat$Var1[which(gene_stat$Freq >5 )])

  for (len in 1:length(test_gene)) {
    DEGMarkers <- test_gene[[len]]
    padj_stat <- c()
    for (clu in 1:ncol(DEGMarkers)) {
      p_stat <- c()
      DEGMarkers_sub <- as.data.frame(DEGMarkers[,clu])
      size_B <- dim(DEGMarkers_sub)[1]
      for (ref in unique(HallMarker1$ref_type1)) {
        HallMarker_sub <- subset(HallMarker1,ref_type1 == ref)
        size_A <- dim(HallMarker_sub)[1]
        size_C <- length(intersect(DEGMarkers_sub[,1],HallMarker_sub$makrer))
        p_value=phyper(size_C-1, size_A, population_size-size_A, size_B, lower.tail=F)
        stat <- data.frame(test=colnames(DEGMarkers)[clu],ref=ref,intersect_num=size_C,p_val=p_value)
        p_stat <- rbind(p_stat,stat)
      }
      padj_stat <- rbind(padj_stat,p_stat)
    }
    set.seed(123)
    padj_stat$p_val_adj <- p.adjust(padj_stat$p_val,method = method)
    #SSEA_list[[paste0(names(test_gene)[len],"_padj_stat")]] <- padj_stat

    padj_stat1 <- padj_stat[,c("test","ref","p_val_adj")]
    colnames(padj_stat1)[3] <- "value"
    padj_stat_mat <- as.matrix(reshape2::dcast(padj_stat1,test~ref))
    rownames(padj_stat_mat) <- padj_stat_mat[,1]
    padj_stat_mat <- padj_stat_mat[,-1]

    min_type <- c()
    for (i in 1:nrow(padj_stat_mat)) {
      pre_padj_stat <- as.data.frame(padj_stat_mat[i,])
      colnames(pre_padj_stat) <- rownames(padj_stat_mat)[i]
      pre_padj_stat$type <- rownames(pre_padj_stat)
      pre_padj_stat[,1] <- as.numeric(pre_padj_stat[,1])
      pre_padj_stat <- pre_padj_stat[order(pre_padj_stat[,1]),]
      min_cols <- data.frame(spot=rownames(padj_stat_mat)[i],p.adj=pre_padj_stat[1,1],predict_type=pre_padj_stat[1,2],top=names(test_gene)[len])
      result_df <- rbind(result_df,min_cols)
    }


  }
  SSEA_list[["predict_df"]] <- as.data.frame(result_df)
  df_result <- c()
  for (i in unique(result_df$spot)) {
    df1 <- subset(result_df,spot==i)
    df_frame <- data.frame(cell=i,predict_spot=names(table(df1$predict_type))[which(table(df1$predict_type) == max(table(df1$predict_type)))])
    if (dim(df_frame)[1] > 1) {
      df_frame_fil <- subset(df1,predict_type %in% unique(df_frame$predict))
      df_frame <- data.frame(cell=i,predict_spot=df_frame_fil$predict_type[which.min(df_frame_fil$p.adj)],padj=df_frame_fil$p.adj[which.min(df_frame_fil$p.adj)])
      if (df_frame$padj < 0.05) {
        df_frame <- df_frame[,c("cell","predict_spot")]
        df_result <- rbind(df_result,df_frame)
      }else{df_frame <- data.frame(cell=i,predict_spot="Others")
      df_result <- rbind(df_result,df_frame)}

    }else{df_frame_fil <- subset(df1,predict_type %in% unique(df_frame$predict))
    if (min(df_frame_fil$p.adj) < 0.05) {
      df_result <- rbind(df_result,df_frame)
    }else{df_frame <- data.frame(cell=i,predict_spot="Others")
    df_result <- rbind(df_result,df_frame)}
    }}
  SSEA_list[["predict_spot"]] <- as.data.frame(df_result)


  ####anno sub type

  result_df_1 <- c()
  df_result_1 <- c()
  for (pre in unique(df_result$predict_spot)) {
    cell <- df_result$cell[which(df_result$predict_spot == pre)]
    HallMarker2 <- HallMarker[which(HallMarker$ref_type1 == pre),c("ref_type2","makrer")]
    gene_stat <- data.frame(table(HallMarker2$ref_type2))
    HallMarker2 <- subset(HallMarker2,ref_type2 %in% gene_stat$Var1[which(gene_stat$Freq >3 )])
    HallMarker2 <- subset(HallMarker2,ref_type2 != "Others")

    if (is.null(dim(HallMarker2))) {
      print("your table of sub type markers is null")
    }else{
      for (len in 1:length(test_gene)) {
        DEGMarkers <- test_gene[[len]]
        if (length(cell) == 1) {
          DEGMarkers <- as.data.frame(DEGMarkers[,cell])
          colnames(DEGMarkers) <- cell
        }else{DEGMarkers <- as.data.frame(DEGMarkers[,cell])}

        padj_stat <- c()
        for (clu in 1:ncol(DEGMarkers)) {
          p_stat <- c()
          DEGMarkers_sub <- as.data.frame(DEGMarkers[,clu])
          size_B <- dim(DEGMarkers_sub)[1]
          for (ref in unique(HallMarker2$ref_type2)) {
            HallMarker_sub <- subset(HallMarker2,ref_type2 == ref)
            size_A <- dim(HallMarker_sub)[1]
            size_C <- length(intersect(DEGMarkers_sub[,1],HallMarker_sub$makrer))
            p_value=phyper(size_C-1, size_A, population_size-size_A, size_B, lower.tail=F)
            stat <- data.frame(test=colnames(DEGMarkers)[clu],ref=ref,intersect_num=size_C,p_val=p_value)
            p_stat <- rbind(p_stat,stat)
          }
          padj_stat <- rbind(padj_stat,p_stat)
        }
        set.seed(123)
        padj_stat$p_val_adj <- p.adjust(padj_stat$p_val,method = method)
        #SSEA_list[[paste0(names(test_gene)[len],"_padj_stat")]] <- padj_stat

        padj_stat1 <- padj_stat[,c("test","ref","p_val_adj")]
        colnames(padj_stat1)[3] <- "value"
        padj_stat_mat <- as.matrix(reshape2::dcast(padj_stat1,test~ref))
        if (dim(padj_stat_mat)[1] > 1) {
          rownames(padj_stat_mat) <- padj_stat_mat[,1]
          padj_stat_mat <- padj_stat_mat[,-1]
        }else{
          rownames(padj_stat_mat) <- padj_stat_mat[,1]
          padj_stat_mat <- t(as.matrix(padj_stat_mat[,-1]))
          rownames(padj_stat_mat) <- cell}


        min_type <- c()
        for (i in 1:nrow(padj_stat_mat)) {
          pre_padj_stat <- as.data.frame(padj_stat_mat[i,])
          colnames(pre_padj_stat) <- rownames(padj_stat_mat)[i]
          pre_padj_stat$type <- rownames(pre_padj_stat)
          pre_padj_stat[,1] <- as.numeric(pre_padj_stat[,1])
          pre_padj_stat <- pre_padj_stat[order(pre_padj_stat[,1]),]
          min_cols <- data.frame(spot=rownames(padj_stat_mat)[i],p.adj=pre_padj_stat[1,1],predict_type=pre_padj_stat[1,2],top=names(test_gene)[len])
          result_df_1 <- rbind(result_df_1,min_cols)
        }
      }
    }
  }
  for (i in unique(result_df_1$spot)) {
    df1 <- subset(result_df_1,spot==i)
    df_frame <- data.frame(cell=i,predict_spot_sub=names(table(df1$predict_type))[which(table(df1$predict_type) == max(table(df1$predict_type)))])
    if (dim(df_frame)[1] > 1) {
      df_frame_fil <- subset(df1,predict_type %in% unique(df_frame$predict))
      df_frame <- data.frame(cell=i,predict_spot_sub=df_frame_fil$predict_type[which.min(df_frame_fil$p.adj)],padj=df_frame_fil$p.adj[which.min(df_frame_fil$p.adj)])
      if (df_frame$padj < 0.05) {
        df_frame <- df_frame[,c("cell","predict_spot_sub")]
        df_result_1 <- rbind(df_result_1,df_frame)
      }else{df_frame <- data.frame(cell=i,predict_spot_sub="Others")
      df_result_1 <- rbind(df_result_1,df_frame)}

    }else{df_frame_fil <- subset(df1,predict_type %in% unique(df_frame$predict))
    if (min(df_frame_fil$p.adj) < 0.05) {
      df_result_1 <- rbind(df_result_1,df_frame)
    }else{df_frame <- data.frame(cell=i,predict_spot_sub="Others")
    df_result_1 <- rbind(df_result_1,df_frame)}
    }}

  SSEA_list[["predict_df_sub"]] <- as.data.frame(result_df_1)
  SSEA_list[["predict_spot_sub"]] <- as.data.frame(df_result_1)
  #SSEA_list[["result_type"]] <- as.data.frame(result_type)
  return(SSEA_list)
}
