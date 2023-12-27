# DrSpace
<!-- badges: start -->
<!-- badges: end -->

DrSpace is ...

## Installation

You can install from Github

``` r
#copykat
if (!requireNamespace("copykat", quietly = TRUE)) { 
    devtools::install_github("navinlabcode/copykat")
}
#CellChat
if (!requireNamespace("CellChat", quietly = TRUE)) { 
    devtools::install_github("jinworks/CellChat")
}
#DrSpace
if (!requireNamespace("CellChat", quietly = TRUE)) { 
    devtools::install_github("QianChwnLyn/DrSpace")
}
```

## Example

示例数据可以在这里下载[数据](https://github.com/QianChwnLyn/DrSpace/tree/main/data)。确保空转数据的images以及clusters的完成。

``` r
library(DrSpace)
Library(Seurat)
load("obj.rda")
Seurat::SpatialDimPlot(obj, pt.size = 1.5,label = TRUE,label.size =3 )

```

利用copykat预测疾病数据

``` r
copy_obj <- Copykat(obj = obj,cancer = "colon cancer",n_PC = 10,genome = "hg20")
Seurat::SpatialDimPlot(copy_obj[[1]], pt.size = 1.5,label = TRUE,label.size =2,group.by = "type")

```

利用ssea对空间转录组数据进行细胞类型富集分析，预测细胞类型。

``` r
num_list <- seq(100,1000,100)
pred_obj <- SSEA(obj_list = copy_obj, num_list, cancer = "colon cancer", population_size = 20000)
anno_obj <- pred_obj[[9]]
Seurat::SpatialDimPlot(anno_obj, pt.size = 1.5,label = TRUE,label.size =2,group.by = "predict_spot")
Seurat::SpatialDimPlot(anno_obj, pt.size = 1.5,label = TRUE,label.size =2,group.by = "predict_spot_sub")
Seurat::SpatialDimPlot(anno_obj, pt.size = 1.5,label = TRUE,label.size =2,group.by = "predict_cluster")

```

利用CellChat对空间转录组数据进行spot-spot communication network 构建。

``` r
ssc_pre <- SSC(anno_obj,json_path = "../data/spatial/scalefactors_json.json")
pathways.show <- “IL6”
CellChat::netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
CellChat::netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5)
CellChat::netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial",edge.width.max = 2, alpha.image = 0.5, vertex.weight = "outgoing", vertex.size.max = 3.5, vertex.label.cex = 3.5)
CellChat::netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial",edge.width.max = 2, alpha.image = 0.5, vertex.weight = "incoming", vertex.size.max = 3.5, vertex.label.cex = 3.5)

```

## References
1.Gao, R Jin et al., Delineating copy number and clonal substructure in human tumors from single-cell transcriptomes. Nat Biotechnol. doi:10.1038/s41587-020-00795-2.

2.Suoqin Jin et al., CellChat for systematic analysis of cell-cell communication from single-cell and spatially resolved transcriptomics, bioRxiv 2023


