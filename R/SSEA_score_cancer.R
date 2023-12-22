#' Spatial genesets enrichment score for cancer spots
#'
#' Predicting cell types in tumor regions from spatial transcriptomic data.
#'
#' @param HallMakrer Gene markers for cancer spots.
#' @param test_gene Top_gene lists.
#' @param population_size Number of genes in the species studied.
#' @param method Statistical methods for p-value correction. c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#'
#' @return A cell type prediction list, containing spot and cluster cell type prediction.
#' @export
#'
#' @examples
#' \dontrun{pred <- SSEA_score_cancer(HallMakrer,test_gene,20000)}
SSEA_score_cancer <- function(HallMakrer,test_gene,population_size,method="BH"){
  colnames(HallMakrer) <- c("ref_type","makrer")
  SSEA_list <- list()
  result_df <- c()
  #result_type <- c()
  for (len in 1:length(test_gene)) {
    DEGMarkers <- test_gene[[len]]
    padj_stat <- c()
    for (clu in 1:ncol(DEGMarkers)) {
      p_stat <- c()
      DEGMarkers_sub <- as.data.frame(DEGMarkers[,clu])
      size_B <- dim(DEGMarkers_sub)[1]
      for (ref in unique(HallMakrer$ref_type)) {
        HallMakrer_sub <- subset(HallMakrer,ref_type == ref)
        size_A <- dim(HallMakrer_sub)[1]
        size_C <- length(intersect(DEGMarkers_sub[,1],HallMakrer_sub$makrer))
        p_value=phyper(size_C-1, size_A, population_size-size_A, size_B, lower.tail=F)
        stat <- data.frame(test=colnames(DEGMarkers)[clu],ref=ref,intersect_num=size_C,p_val=p_value)
        p_stat <- rbind(p_stat,stat)
      }
      padj_stat <- rbind(padj_stat,p_stat)
    }
    padj_stat$p_val_adj <- p.adjust(padj_stat$p_val,method = method)
    #SSEA_list[[paste0(names(test_gene)[len],"_padj_stat")]] <- padj_stat

    padj_stat1 <- padj_stat[,c("test","ref","p_val_adj")]
    colnames(padj_stat1)[3] <- "value"
    padj_stat_mat <- as.matrix(reshape2::dcast(padj_stat1,test~ref))
    rownames(padj_stat_mat) <- padj_stat_mat[,1]
    padj_stat_mat <- padj_stat_mat[,-1]

    min_type <- c()
    for (i in 1:nrow(padj_stat_mat)) {
      #padj_stat_mat[i,] <- as.numeric(padj_stat_mat[i,])
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
  #SSEA_list[["result_type"]] <- as.data.frame(result_type)
  return(SSEA_list)
}
