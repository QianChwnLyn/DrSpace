#' Spatial Gene Sets Enrichment Analysis
#'
#' Separate annotation for tumor spots and normal spots processed data. Integrate tumor and normal control predictions into original object.
#'
#' @param obj_list A seurat object list. Includes whole spatial transcriptome data and data split by disease and normal controls.
#' @param num_list A numeric list. e.g. seq(100,1000,100).
#' @param cancer Cancer type.
#' @param population_size Number of genes in the species studied.
#' @param method Statistical methods for p-value correction. c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#'
#' @return A list, containing spot and cluster cell type prediction and predicted object.
#' @export
#'
#' @examples
#' \dontrun{pred <- SSEA(obj_list,num_list,"colon cancer",20000)}
SSEA <- function(obj_list,num_list,cancer=cancer,population_size,method="BH"){

  ssea_result <- list()
  ###cancer
  cancer_obj <- obj_list[[cancer]]
  mat <- cancer_obj@assays$Spatial@counts
  meta <- data.frame(cell=rownames(cancer_obj@meta.data),cluster=cancer_obj$seurat_clusters)
  ####calculate top gene
  test_gene <- top_gene(mat,num_list)
  wb <- openxlsx::createWorkbook()
  for (i in seq_along(test_gene)) {
    sheet_name <- names(test_gene)[i]
    openxlsx::addWorksheet(wb,sheetName = sheet_name)
    openxlsx::writeData(wb,sheet = sheet_name,x=test_gene[[i]])
  }
  openxlsx::saveWorkbook(wb,paste0(cancer,"_topgene.xlsx"),overwrite = TRUE)

    HallMakrer <- cancer_list
    pred <- SSEA_score_cancer(HallMakrer,test_gene,population_size,method=method)
    cancer_df <- pred[[2]]
    meta <- merge(meta,cancer_df,"cell")
    predict_table <- data.frame(table(meta$cluster,meta$predict_spot))
    predict_table <- reshape2::dcast(predict_table,Var1~Var2)
    rownames(predict_table) <- predict_table[,1]
    predict_table <- predict_table[,-1]

    ###calculate ratio
    predict_cluster <- c()
    predict_table_per <- t(apply(predict_table, 1, function(x) if (sum(x) != 0) {round(x/sum(x)*100)}))
    for (i in 1:nrow(predict_table_per)) {
      predict_table_per_fil <- as.data.frame(predict_table_per[i,])
      colnames(predict_table_per_fil) <- "percentage"
      predict_table_per_fil <- subset(predict_table_per_fil,percentage > 5)
      predict_table_per_fil$type <- paste0(predict_table_per_fil$percentage,"%",rownames(predict_table_per_fil))
      predict_table_per_fil <- predict_table_per_fil[order(predict_table_per_fil$percentage,decreasing = T),]
      predict_table_per_fil_stat <- data.frame(cluster=as.numeric(rownames(predict_table_per)[i]), predict_cluster = paste(predict_table_per_fil$type,collapse = "+"))
      predict_cluster <- rbind(predict_cluster,predict_table_per_fil_stat)
    }

    meta$predict_cluster <- as.numeric(as.character(meta$cluster))
    for (i in unique(meta$predict_cluster)) {
      meta$predict_cluster[which(meta$predict_cluster == i)] <- predict_cluster$predict_cluster[which(predict_cluster$cluster == i)]
    }
    pred[["predict_cancer_cluster"]] <- meta
    ssea_result[1:3] <- pred
    ###normal
    normal <- obj_list[["normal"]]
    mat <- normal@assays$Spatial@counts
    meta <- data.frame(cell=rownames(normal@meta.data),cluster=normal$seurat_clusters)
    ####calculate top gene
    test_gene <- top_gene(mat,num_list)
    wb <- openxlsx::createWorkbook()
    for (i in seq_along(test_gene)) {
      sheet_name <- names(test_gene)[i]
      openxlsx::addWorksheet(wb,sheetName = sheet_name)
      openxlsx::writeData(wb,sheet = sheet_name,x=test_gene[[i]])
    }
    openxlsx::saveWorkbook(wb,paste0("normal","_topgene.xlsx"),overwrite = TRUE)
    method=method

    HallMakrer <- immune_list
    pred1 <- SSEA_score_normal(HallMakrer,test_gene,population_size,method=method)
    normal_df <- pred1[[2]]
    normal_df1 <- pred1[[4]]
    normal_df <- merge(normal_df,normal_df1,"cell")

    meta <- merge(meta,normal_df,"cell")
    meta$predict_spot_sub <- paste0(meta$predict_spot,": ",meta$predict_spot_sub)
    predict_table <- data.frame(table(meta$cluster,meta$predict_spot_sub))
    predict_table <- reshape2::dcast(predict_table,Var1~Var2)
    rownames(predict_table) <- predict_table[,1]
    predict_table <- predict_table[,-1]

    ###calculate ratio
    predict_cluster <- c()
    predict_table_per <- t(apply(predict_table, 1, function(x) if (sum(x) != 0) {round(x/sum(x)*100)}))
    for (i in 1:nrow(predict_table_per)) {
      predict_table_per_fil <- as.data.frame(predict_table_per[i,])
      colnames(predict_table_per_fil) <- "percentage"
      predict_table_per_fil <- subset(predict_table_per_fil,percentage > 5)
      predict_table_per_fil$type <- paste0(predict_table_per_fil$percentage,"%",rownames(predict_table_per_fil))
      predict_table_per_fil <- predict_table_per_fil[order(predict_table_per_fil$percentage,decreasing = T),]
      predict_table_per_fil_stat <- data.frame(cluster=as.numeric(rownames(predict_table_per)[i]), predict_cluster = paste(predict_table_per_fil$type,collapse = "+"))
      predict_cluster <- rbind(predict_cluster,predict_table_per_fil_stat)
    }

    meta$predict_cluster <- as.numeric(as.character(meta$cluster))
    for (i in unique(meta$predict_cluster)) {
      meta$predict_cluster[which(meta$predict_cluster == i)] <- predict_cluster$predict_cluster[which(predict_cluster$cluster == i)]
    }
    pred1[["predict_normal_cluster"]] <- meta
    ssea_result[4:8] <- pred1
    names(ssea_result) <- c(names(pred),names(pred1))
    ###integrate normal+cancer
    df_normal <- pred1[[5]]
    df_cancer <- pred[[3]]
    df_cancer$predict_spot_sub <- df_cancer$predict_spot
    df_cancer <- df_cancer[c(1:3,5,4)]
    head(df_normal)
    head(df_cancer)
    df_result <- rbind(df_normal,df_cancer)

    obj <- obj_list[["obj"]]
    obj$cell <- rownames(obj@meta.data)
    rownames(df_result) <- df_result$cell
    obj@meta.data <- merge(obj@meta.data,df_result[,c(1,3:5)])
    rownames(obj@meta.data) <- obj$cell
    ssea_result[["predict_obj"]] <- obj


  return(ssea_result)
}
