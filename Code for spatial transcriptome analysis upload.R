library(Seurat)
library(future)
library(remotes)
library(ggplot2)
library(patchwork)
library(dplyr)
library(progeny)
library(tidyr)
library(readr)
library(pheatmap)
library(tibble)
library(viridis)

source('Bmk_Space_mapping.R')
source("CreateBmkObject.R")

i= 3 #选择spots尺寸等级
object_level3 <- CreateS1000Object(
  matrix_path=paste0('D:/XXX/XXX/矩阵文件路径',i,'_heAuto'),#矩阵文件目录
  png_path="HCCNC1_dyeing.png", #png格式图片路径
  spot_radius =1, #点的半径，可以不指定，会自动计算
  min.cells = 10, #一个基因至少在n个细胞中表达才被保留，可自行调整，默认值5
  min.features = 0 #一个细胞至少有n个基因才被保留，可自行调整，默认值100
)

#SCT标准化
object_level3 <- SCTransform(object_level3, assay = "Spatial", verbose = FALSE)

#Gene expression visualization
SpatialFeaturePlot(object_level3, features = "FOS",pt.size.factor = 0.009,image.alpha = 0.5,alpha = c(0.5,1),min.cutoff = 'q01',
                   max.cutoff = 'q99') + theme(legend.position = "right") 

plot <- SpatialFeaturePlot(object_level3_subset, features = c("Gja5"),pt.size.factor = 0.003) + theme(legend.text = element_text(size = 0),
                                                                                                            legend.title = element_text(size = 20), legend.key.size = unit(1, "cm"))
jpeg(filename = "spatial_vignette_ttr.jpg", height = 700, width = 1200, quality = 50)
print(plot)
dev.off()

#Dimensionality reduction, clustering, and visualization
object_level3 <- RunPCA(object_level3, assay = "SCT", verbose = FALSE)
object_level3 <- FindNeighbors(object_level3, reduction = "pca", dims = 1:30)
object_level3 <- FindClusters(object_level3, verbose = T,resolution=4)
object_level3 <- RunUMAP(object_level3, reduction = "pca", dims = 1:30)
DimPlot(object_level3, reduction = "umap", label = TRUE)
SpatialDimPlot(object_level3, label = TRUE,pt.size.factor = 0.0035, label.size = 3)

SpatialDimPlot(object_level3, cells.highlight = CellsByIdentities(object = object_level3, idents = c(0:12)), facet.highlight = TRUE, ncol = 4,pt.size.factor = 0.0035)

#整合单细胞
library(dplyr)
scRNA_DATA <- SCTransform(scRNA_DATA, ncells = 3000, verbose = T) %>%
  RunPCA(verbose = T) %>%
  RunUMAP(dims = 1:30)

object_level3 <- SCTransform(object_level3, assay = "Spatial", verbose = T) %>%
  RunPCA(verbose = T)
# the annotation is stored in the 'subclass' column of object metadata
DimPlot(scRNA_DATA, group.by = "seurat_7clusters", label = TRUE)

anchors <- FindTransferAnchors(reference = scRNA_DATA, query = object_level3, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = scRNA_DATA$seurat_7clusters, prediction.assay = TRUE,
                                  weight.reduction = object_level3[["pca"]], dims = 1:30)
object_level3[["predictions"]] <- predictions.assay
#映射位置
DefaultAssay(object_level3) <- "predictions"
SpatialFeaturePlot(object_level3, features = c(""), pt.size.factor = 0.012,
                   ncol = 1, crop = TRUE,image.alpha = 0.5,min.cutoff = 'q02', max.cutoff = 'q98')+scale_fill_gradient2(low="#372f76",mid="#ae465b",high="yellow",midpoint=0.4)+ theme(legend.position = "right")

#基于差异基因映射位置
object_level3 <- FindSpatiallyVariableFeatures(object_level3, assay = "predictions", selection.method = "moransi", r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(object_level3, selection.method = "moransi"), 1)
SpatialPlot(object = object_level3, features = top.clusters, ncol = 3, pt.size.factor = 0.007, crop = TRUE,image.alpha = 0.5,min.cutoff = 'q02',
            max.cutoff = 'q98')+scale_fill_continuous(low="#eeeeee",high="blue")
#查看预测类型结果
object_level3$predicted.id <- GetTransferPredictions(object_level3)
Idents(object_level3) <- "predicted.id"
table(Idents(object_level3))

#progeny信号通路
object_level3<-progeny(object_level3,scale = TRUE,organism = "Human",top = 500,assay_name = "Spatial",return_assay = T)
object_level3 <- Seurat::ScaleData(object_level3, assay = "progeny")
object_level3@meta.data<-cbind(object_level3@meta.data,t(object_level3@assays$progeny@scale.data))
#输出图像
listfori<-unique(rownames(object_level3@assays$progeny@data))
listfori<-c("p53","PI3K","TGFb")
for (i in listfori) {
  plot <- SpatialPlot(object = object_level3, features = i, ncol = 1, pt.size.factor = 0.006,
                      crop = TRUE,image.alpha = 1,min.cutoff = 'q02', max.cutoff = 'q98')+scale_fill_gradient2(low="#372f76",mid="#21908CFF",high="yellow",midpoint = 1.5)+ theme(legend.position = "right")
  jpeg(filename = paste("object_level3",i,".jpeg"), height = 4800, width = 6400, quality = 100,res=600)  
  print(plot)
  dev.off()
}

listfori<-unique(rownames(object_level3@assays$progeny@data))
listfori<-c("p53","PI3K","TGFb")
for (i in listfori) {
  SpatialFeaturePlot(object = object_level3, features = "VWF", ncol = 1, pt.size.factor = 0.012,
                     crop = TRUE,image.alpha = 0.5,min.cutoff = 'q02', max.cutoff = 'q98',alpha = c(0.1,1))+scale_fill_gradient(low="#372f76",high="green")+ theme(legend.position = "right")
}