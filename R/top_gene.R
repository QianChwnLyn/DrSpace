#' Top n genes
#'
#' According to the provided numeric list, stat the genes whose expression level is in the top n for each spot.
#'
#' @param mat A gene expression matrix. The row name is the gene name and the column name is spot.
#' @param n A numeric list. e.g. seq(100,1000,100).
#'
#' @return The default method returns a data frame list.
#' @export
#'
#' @examples
#' \dontrun{top_gene(mat,num_list)}
top_gene <- function(mat,n){
  gene_top_list <- list()
  for (num in n) {
    matr <- as.data.frame(mat[,1])
    colnames(matr) <- "cell"
    matr$gene <- rownames(matr)
    matr <- as.data.frame(matr[order(matr$cell,decreasing = T),])
    gene_table <- as.data.frame(matr$gene[1:num])
    colnames(gene_table) <- colnames(mat)[1]
    for (spot in 2:ncol(mat)) {
      matr <- as.data.frame(mat[,spot])
      colnames(matr) <- "cell"
      matr$gene <- rownames(matr)
      matr <- as.data.frame(matr[order(matr$cell,decreasing = T),])
      gene <- as.data.frame(matr$gene[1:num])
      colnames(gene) <- colnames(mat)[spot]
      gene_table <- cbind(gene_table,gene)
    }
    gene_top_list[[paste0("top_",num)]] <- gene_table
    #write.table(gene_table,paste0("top_",n,"_gene.txt"),quote = F,sep = "\t",row.names = F)
  }
  return(gene_top_list)
}
