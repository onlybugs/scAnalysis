# 22-12-10
# kli
# ... ST rna

# install.packages('hdf5r')

library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(hdf5r)

#########################################################################
# 读数据 可以先读矩阵
expr <- "./visual/Visium_FFPE_Human_Breast_Cancer_filtered_feature_bc_matrix.h5"
expr.data <- Read10X_h5(filename = expr)
mydata <- CreateSeuratObject(counts = expr.data,project = 'HBC',assay = 'Spatial')
mydata$slice <- 1
mydata$region <- 'HBC'
# 读图
img <- Read10X_Image("./visual/Visium_FFPE_Human_Breast_Cancer_spatial/spatial")
DefaultAssay(img) <- 'Spatial'
img <- img[colnames(x = mydata)]
mydata[['image']] <- img

# 直接读的方法 其实差不多 推荐上面的 可操作性更强
hbc <- Load10X_Spatial(data.dir = "./visual/Visium_FFPE_Human_Breast_Cancer_spatial",
                       filename = "Visium_FFPE_Human_Breast_Cancer_filtered_feature_bc_matrix.h5",
                       assay = 'Spatial',
                       slice = 'spatial',
                       filter.matrix = T)

#########################################################################
p1 <- VlnPlot(mydata,features = "nCount_Spatial",pt.size = 0.1) + NoLegend()
p2 <- SpatialFeaturePlot(mydata,features = "nCount_Spatial") + theme(legend.position = 'right')
plot_grid(p1,p2)

mydata[['percent.mt']] <- PercentageFeatureSet(mydata,pattern = "^mt[-]")
p1 <- VlnPlot(mydata,features = "percent.mt",pt.size = 0.1) + NoLegend()
p2 <- SpatialFeaturePlot(mydata,features = "percent.mt") + theme(legend.position = 'right')
plot_grid(p1,p2)

mydata <- subset(mydata,
               subset = nFeature_Spatial>200 & nFeature_Spatial < 7500 &
                 percent.mt < 20 & nCount_Spatial < 6000)
# 标准化
mydata <- SCTransform(mydata,assay = "Spatial")
# 特定基因表达
SpatialFeaturePlot(mydata,features = c("TP53","EGFR"))

# PCA UMAP
mydata <- RunPCA(mydata,assay = "SCT")
mydata <- FindNeighbors(object = mydata,dims = 1:8)
mydata <- FindClusters(mydata,resolution = 0.5)
mydata <- RunUMAP(mydata,reduction = 'pca',dims = 1:8)

p1 <- DimPlot(mydata,reduction = 'umap',label = T)
p2 <- SpatialDimPlot(mydata,label = T,label.size = 3)
plot_grid(p1,p2)


# 余下分析差不多



