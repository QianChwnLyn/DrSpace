#' Spatial transcriptomic cancer spot prediction
#'
#' Use copykat to predict cancer spots and split normal spots and cancer spots.
#'
#' @param obj Spatial transcriptome object.
#' @param cancer Cancer type.
#' @param n_PC Dimensions of reduction to use as input.
#' @param id.type Gene id type: Symbol or Ensemble.
#' @param cell.line If the data are from pure cell line,put "yes"; if cell line data are a mixture of tumor and normal cells, still put "no".
#' @param ngene.chr Minimal number of genes per chromosome for cell filtering.
#' @param LOW.DR Minimal population fractions of genes for smoothing.
#' @param UP.DR Minimal population fractions of genes for segmentation.
#' @param win.size Minimal window sizes for segmentation.
#' @param norm.cell.names A vector of normal cell names.
#' @param KS.cut Segmentation parameters, input 0 to 1; larger looser criteria.
#' @param sam.name Sample name.
#' @param distance Distance methods include euclidean, and correlation converted distance include pearson and spearman.
#' @param output.seg TRUE or FALSE, output seg file for IGV visualization.
#' @param plot.genes TRUE or FALSE, output heatmap of CNVs with genename labels.
#' @param genome hg20 or mm10, current version only work for human or mouse genes.
#' @param n.cores Number of cores for parallel computing.
#'
#' @return A list. containing predicted object, cancer object and normal object.
#' @export
#'
#' @examples
#' \dontrun{test <- Copykat(obj,cancer = "colon cancer",n_PC = 10,id.type="S", ngene.chr=5, win.size=25, KS.cut=0.1, sam.name="GSM7058756_C1", distance="euclidean", norm.cell.names="")}
Copykat <- function(obj,cancer="cancer",n_PC=10,id.type = "S",cell.line = "no",ngene.chr = 5,LOW.DR = 0.05,UP.DR = 0.1,
                    win.size = 25,norm.cell.names = "",KS.cut = 0.1,sam.name = "",distance = "euclidean",output.seg = "FALSE",plot.genes = "TRUE",genome = "hg20",n.cores = 1){
  obj_list <- list()
  if (!requireNamespace("copykat",quietly = TRUE)) {
    devtools::install_github("navinlabcode/copykat")
  }
  library(copykat)
  exp.rawdata <- as.matrix(obj@assays$Spatial@counts)
  copykat.test <- copykat::copykat(rawmat=exp.rawdata, id.type=id.type, cell.line = cell.line, ngene.chr = ngene.chr, LOW.DR = LOW.DR, UP.DR = UP.DR, win.size =win.size,
                          norm.cell.names =norm.cell.names, KS.cut = KS.cut, sam.name = sam.name, distance = distance, output.seg = output.seg, plot.genes = plot.genes, genome = genome, n.cores = n.cores)
  pred.test <- data.frame(copykat.test$prediction)
  saveRDS(copykat.test,"copykat.test.rds")
  colnames(pred.test)[1] <- "Cell"
  obj$Cell <- rownames(obj@meta.data)
  obj@meta.data <- dplyr::left_join(obj@meta.data,pred.test,"Cell")
  rownames(obj@meta.data) <- obj$Cell
  obj$type <- "Normal"
  obj$type[which(obj$copykat.pred == "aneuploid")] <- cancer
  obj$type <- factor(obj$type,levels = c(cancer,"Normal"))
  Idents(obj) <- obj$type
  obj_list[["obj"]] <- obj
  ###subset normal
  normal <- subset(obj,type=="Normal")
  pdf("NormalType_subset.pdf",width = 12,height = 6)
  p1 <- Seurat::SpatialDimPlot(normal, pt.size = 1.5,label = TRUE,label.size =3 ,group.by = "type")
  p2 <- Seurat::SpatialFeaturePlot(normal,features = "EPCAM",pt.size = 1.5)
  print(patchwork::wrap_plots(p1,p2))
  dev.off()
  normal<-Seurat::RunPCA(normal,verbose = FALSE)
  normal<-Seurat::FindNeighbors(normal,dims = 1:n_PC)
  normal<-Seurat::FindClusters(normal,resolution = 2)
  obj_list[["normal"]] <- normal
  ###subset cancer
  obj1 <- subset(obj,type==cancer)
  pdf("CancerType_subset.pdf",width = 12,height = 6)
  p1 <- Seurat::SpatialDimPlot(obj1, pt.size = 2,label = TRUE,label.size =3 ,group.by = "type")
  p2 <- Seurat::SpatialFeaturePlot(obj1,features = "EPCAM",pt.size = 2)
  print(patchwork::wrap_plots(p1,p2))
  dev.off()
  obj1<-Seurat::RunPCA(obj1,verbose = FALSE)
  obj1<-Seurat::FindNeighbors(obj1,dims = 1:n_PC)
  obj1<-Seurat::FindClusters(obj1,resolution = 2)
  obj_list[[cancer]] <- obj1

  return(obj_list)

}
