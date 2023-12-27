# DrSpace: Disease-rating based on spatially resolved cancer micro-environment
<!-- badges: start -->
<!-- badges: end -->

#### 1. Description

**Aim to solve**: 

* Annotate tumor micro-environment:
  * Achieve distinct annotation of malignant cells and immune cells in the tumor microenvironment.
  * Enable precise annotation of immune cell subtypes.

* On the basis of annotating tumor micro-environment, we employ `CellChat` to construct a tumor micro-ecological network.

#### 2. Installation

Installing `Copykat` from GitHub.

```R
if (!requireNamespace("copykat", quietly = TRUE)) { 
    devtools::install_github("navinlabcode/copykat")
}
```

Installing `CellChat ` from GitHub.

```r
if (!requireNamespace("CellChat", quietly = TRUE)) { 
    devtools::install_github("jinworks/CellChat")
}
```

Installing `DrSpace` from GitHub.

```r
if (!requireNamespace("CellChat", quietly = TRUE)) { 
    devtools::install_github("QianChwnLyn/DrSpace")
}
```

#### 3. Usage

Example data can be downloaded [here](https://github.com/QianChwnLyn/DrSpace/tree/main/data). Make sure to complete the images and clusters for the empty revolving data.

```R
library(DrSpace)
library(Seurat)
load("obj.rda")
Seurat::SpatialDimPlot(obj, pt.size = 1.5,label = TRUE,label.size =3 )
```

Predict disease data using `Copykat`.

```r
copy_obj <- Copykat(obj = obj,cancer = "colon cancer",n_PC = 10,genome = "hg20")
Seurat::SpatialDimPlot(copy_obj[[1]], pt.size = 1.5,label = TRUE,label.size =2,group.by = "type")
```

Perform cell type enrichment analysis and predicted cell types on spatial transcriptomic data using `SSEA`.

```r
num_list <- seq(100,1000,100)
pred_obj <- SSEA(obj_list = copy_obj, num_list, cancer = "colon cancer", population_size = 20000)
anno_obj <- pred_obj[[9]]
Seurat::SpatialDimPlot(anno_obj, pt.size = 1.5,label = TRUE,label.size =2,group.by = "predict_spot")
Seurat::SpatialDimPlot(anno_obj, pt.size = 1.5,label = TRUE,label.size =2,group.by = "predict_spot_sub")
Seurat::SpatialDimPlot(anno_obj, pt.size = 1.5,label = TRUE,label.size =2,group.by = "predict_cluster")
```

Construct Spot-Spot Communication Network on spatial transcriptomic data using `CellChat`.

```R
ssc_pre <- SSC(anno_obj,json_path = "../data/spatial/scalefactors_json.json")
pathways.show <- “IL6”
CellChat::netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
CellChat::netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5)
CellChat::netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial",edge.width.max = 2, alpha.image = 0.5, vertex.weight = "outgoing", vertex.size.max = 3.5, vertex.label.cex = 3.5)
CellChat::netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial",edge.width.max = 2, alpha.image = 0.5, vertex.weight = "incoming", vertex.size.max = 3.5, vertex.label.cex = 3.5)
```

#### Reference

1. Gao, R Jin et al., Delineating copy number and clonal substructure in human tumors from single-cell transcriptomes. Nat Biotechnol. doi:10.1038/s41587-020-00795-2.

2. Suoqin Jin et al., CellChat for systematic analysis of cell-cell communication from single-cell and spatially resolved transcriptomics, bioRxiv 2023


