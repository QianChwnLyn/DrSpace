#' Spatial spot-spot communication
#'
#' Using cellchat to calculate spatial transcriptome spot-spot communication network.
#'
#' @param obj Object containing the cell type.
#' @param json_path The path of json.
#' @param spot.size The theoretical spot size (um) in 10X Visium
#'
#' @return The result of spot-spot communication.
#' @export
#'
#' @examples
#' \dontrun{ssc_pre <- SSC(obj,json_path = "../data/spatial/scalefactors_json.json")}
SSC <- function(obj,json_path,spot.size=65){
  Seurat::Idents(obj) <- obj$predict_cluster
  #矩阵信息
  color.use <- CellChat::scPalette(nlevels(obj))
  names(color.use) <- levels(obj)
  data.input = Seurat::GetAssayData(obj, slot = "data", assay = "Spatial")
  #meta信息
  meta = data.frame(labels = Seurat::Idents(obj), slices = "slice1", row.names = names(Seurat::Idents(obj))) # manually create a dataframe consisting of the cell labels
  unique(meta$labels)
  meta$slices <- factor(meta$slices)
  # 空间图像信息
  spatial.locs = Seurat::GetTissueCoordinates(obj, scale = NULL, cols = c("imagerow", "imagecol"))
  # Scale factors and spot diameters 信息
  scalefactors = jsonlite::fromJSON(txt = file.path(json_path))
  #spot.size = 65 # the theoretical spot size (um) in 10X Visium
  conversion.factor = spot.size/scalefactors$spot_diameter_fullres
  spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)
  cellchat <- CellChat::createCellChat(object = data.input, meta = meta, group.by = "labels",
                             datatype = "spatial", coordinates = spatial.locs, spatial.factors = spatial.factors)

  CellChatDB <- CellChat::CellChatDB.human
  CellChatDB.use <- CellChatDB
  cellchat@DB <- CellChatDB.use

  cellchat <- CellChat::subsetData(cellchat)
  future::plan("multisession", workers = 4)
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)


  ptm = Sys.time()
  cellchat <- CellChat::computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                                distance.use = TRUE, interaction.range = 250, scale.distance = 0.01,
                                contact.dependent = TRUE, contact.range = 100)
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = 10)
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)

  execution.time = Sys.time() - ptm
  print(as.numeric(execution.time, units = "secs"))

  ptm = Sys.time()

  groupSize <- as.numeric(table(cellchat@idents))
  pdf("cellchat_net_number.pdf",width = 16,height = 16)
  p1 <- CellChat::netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")
  print(p1)
  dev.off()
  pdf("cellchat_net_weight.pdf",width = 16,height = 16)
  p1 <- CellChat::netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Blues")
  print(p1)
  dev.off()
  pathways.show <- cellchat@netP$pathways
  # Circle plot
  cellchat <- CellChat::netAnalysis_computeCentrality(cellchat, slot.name = "netP")

  for (i in pathways.show) {
    pdf(paste0(i,"_circle.pdf"),width = 6,height = 6)
    p1 <- CellChat::netVisual_aggregate(cellchat, signaling = i, layout = "circle")
    print(p1)
    dev.off()

    pdf(paste0(i,"_spatial.pdf"),width = 10,height = 10)
    p1 <- CellChat::netVisual_aggregate(cellchat, signaling = i, layout = "spatial", edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5)
    print(p1)
    dev.off()


    pdf(paste0(i,"_outgoing_network.pdf"),width = 10,height = 10)
    p1 <- CellChat::netVisual_aggregate(cellchat, signaling = i, layout = "spatial",
                              edge.width.max = 2, alpha.image = 0.5, vertex.weight = "outgoing", vertex.size.max = 3.5, vertex.label.cex = 3.5)
    print(p1)
    dev.off()
    pdf(paste0(i,"_incoming_network.pdf"),width = 10,height = 10)
    p2 <- CellChat::netVisual_aggregate(cellchat, signaling = i, layout = "spatial",
                              edge.width.max = 2, alpha.image = 0.5, vertex.weight = "incoming", vertex.size.max = 3.5, vertex.label.cex = 3.5)
    print(p2)
    dev.off()
  }
  return(cellchat)
}
