# 22-12-8 - 23-1-1
# kli
# scRNA analysis

# install.packages("Seurat")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("GSVA")
# BiocManager::install("GSEABase")
# BiocManager::install("limma")
# BiocManager::install("SingleR")
# BiocManager::install("celldex")
# BiocManager::install("monocle")

library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(celldex)
library(SingleR)
library(monocle)

options(scipen = 200)

###########################1 数据读取和预处理############################
# Read Data 1 from 3 files
data <- Read10X(data.dir = "./pbmc3k_filtered_gene_bc_matrices")
# data[c('TP53',"CD3D"),1:10]
# 创建S对象
pbmc <- CreateSeuratObject(counts = data,project = "pbmc",min.cells = 3,min.features = 1)
head(pbmc@meta.data)
# 线粒体基因
pbmc[['percent.mt']] <- PercentageFeatureSet(pbmc,pattern = "^MT-")
# 提琴图 这个加了meta信息就要使用meta 不要用原始的那个了
VlnPlot(pbmc,
        features = c("nFeature_RNA","nCount_RNA","percent.mt"))
# 绘制任意两个向量的相关性图 所有的值都存储于meta对象中
FeatureScatter(pbmc,feature1 = 'nCount_RNA',
               feature2 = "nFeature_RNA")
# 数据过滤(质控)
pbmc <- subset(pbmc,
               subset = nFeature_RNA>200 & nFeature_RNA < 3000 & percent.mt < 20)

# 标准化
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", 
                         scale.factor = 10000)
# 表达矩阵在这里
pbmc[["RNA"]]@data


# 第二种读取方法 从一个csv文件中进行读取 
rt=read.table("lung/all.csv.gz", header=T, sep=",")
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),
            nrow=nrow(exp),
            dimnames=dimnames)
data=avereps(data)
rm(exp,rt)
lung=CreateSeuratObject(counts = data,
                        project = "lung", 
                        min.cells=3, min.features=50, names.delim = "_")

# 质控
lung[['percent.mt']] <- PercentageFeatureSet(lung,pattern = "^MT-")
# 提琴图 这个加了meta信息就要使用meta 不要用原始的那个了
VlnPlot(lung,
        features = c("nFeature_RNA","nCount_RNA","percent.mt"))
# 数据过滤(质控)
lung <- subset(lung,
               subset = nFeature_RNA>200 & nFeature_RNA < 3000 & percent.mt < 20)

# 标准化
lung <- NormalizeData(lung, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
# 表达矩阵在这里
lung[["RNA"]]@data


###########################2 添加meta数据############################
# 这里去SRA里去找
lun_meta <- read.table("lung/SraRunTable.txt",sep = ",",header = T)
# table(lun_meta$Group) meta数据需要格式
# row name是输入矩阵的列名 即细胞是对应的
# 其他每一列都是一组meta数据
new.lun_meta <- lun_meta[,c("Age",'Group')]
new.lun_meta$sample <- colnames(data)
sce_meta <- data.frame(group = new.lun_meta$Group,
                       row.names = new.lun_meta$sample)
lung=CreateSeuratObject(counts = data,
                        project = "lung", 
                        meta.data = sce_meta,
                        min.cells=3, min.features=50, names.delim = "_")
# head(lung@meta.data)
rm(data,dimnames,lun_meta,sce_meta,new.lun_meta)
# VlnPlot(lung,group.by = "group",
#         features = c("nFeature_RNA","nCount_RNA","percent.mt"))

###########################3 高变基因 降维 聚类############################
# 高变基因
#提取细胞间变异系数较大的基因
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 1500)
#输出特征方差图 并 加入高变基因
top10 <- head(x = VariableFeatures(object = pbmc), 10)
p1 = VariableFeaturePlot(object = pbmc)
LabelPoints(plot = p1, points = top10)
# CombinePlots(plots = list(plot1, plot2))

# 降维
pbmc <- ScaleData(pbmc,features = rownames(pbmc))
pbmc <- RunPCA(pbmc,features = VariableFeatures(object = pbmc))
# 这玩意意义不大
VizDimLoadings(object = pbmc, dims = 1:4, 
               reduction = "pca",nfeatures = 20)
# pc1和pc2的降维图
DimPlot(pbmc,reduction = 'pca')
# 不错的热图
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE,nfeatures = 20)
# 选定PC 都是7-8个合适 pbmc那就八个
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(object = pbmc, dims = 1:20)
ElbowPlot(pbmc,reduction = 'pca')

# 聚类分型了 TSNE 这里就要选降维数了 用8 
pbmc <- FindNeighbors(object = pbmc,dims = 1:7)
pbmc <- FindClusters(pbmc,resolution = 0.5)
# Idents(pbmc) 这个可以查看细胞的分组
pbmc <- RunTSNE(pbmc,dims = 1:7)
TSNEPlot(object = pbmc, pt.size = 2, label = TRUE)

###########################4 差异基因鉴定 marker基因############################
##查找每个聚类的差异基因
pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = 1)
sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>1 &
                            as.numeric(as.vector(pbmc.markers$p_val_adj))<0.05),]
# 如果先进行了注释 那这里就不只是cluster了 还会有其他的
# write.table(sig.markers,file="03.clusterMarkers.txt",sep="\t",row.names=F,quote=F)

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
#绘制marker在各个cluster的热图
DoHeatmap(object = pbmc, features = top10$gene) + NoLegend()

#绘制某些marker基因的小提琴图
VlnPlot(object = pbmc, features = row.names(sig.markers)[1:2])

#需要展示的基因，可以修改
showGenes = top10$gene[1:10]
#绘制marker在各个cluster的散点图
FeaturePlot(object = pbmc, features = showGenes, cols = c("green", "red"))

#绘制marker在各个cluster的气泡图
cluster10Marker=showGenes
DotPlot(object = pbmc, features = cluster10Marker)

# 这里其实是有问题的 最好这个类别不使用聚类的 使用meta数据的更有意义


###########################5 singleR注释############################
counts<-pbmc@assays$RNA@counts
clusters<-pbmc@meta.data$seurat_clusters
ann=pbmc@meta.data$orig.ident
#ref=get(load("ref_Human_all.RData"))
ref=celldex::HumanPrimaryCellAtlasData()
singler=SingleR(test=counts, ref =ref,
                labels=ref$label.main, clusters = clusters)
clusterAnn=as.data.frame(singler)
clusterAnn=cbind(id=row.names(clusterAnn), clusterAnn)
clusterAnn=clusterAnn[,c("id", "labels")]
# write.table(clusterAnn,file="04.clusterAnn.txt",quote=F,sep="\t", row.names=F)
singler2=SingleR(test=counts, ref =ref, 
                 labels=ref$label.main)
cellAnn=as.data.frame(singler2)
cellAnn=cbind(id=row.names(cellAnn), cellAnn)
cellAnn=cellAnn[,c("id", "labels")]
# write.table(cellAnn, file="04.cellAnn.txt", quote=F, sep="\t", row.names=F)

#cluster注释后的可视化
newLabels=singler$labels
names(newLabels)=levels(pbmc)
pbmc=RenameIdents(pbmc, newLabels)
# pdf(file="04.TSNE.pdf",width=6.5,height=6)
TSNEPlot(object = pbmc, pt.size = 2, label = TRUE)    #TSNE可视化
# dev.off()

##cluster注释后的差异分析
pbmc.markers=FindAllMarkers(object = pbmc,
                            only.pos = FALSE,
                            min.pct = 0.25,
                            logfc.threshold = 1)
sig.cellMarkers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>1 & as.numeric(as.vector(pbmc.markers$p_val_adj))<0.05),]
# write.table(sig.cellMarkers,file="04.cellMarkers.txt",sep="\t",row.names=F,quote=F)
# top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
# #绘制marker在各个cluster的热图
# DoHeatmap(object = pbmc, features = top10$gene) + NoLegend()


###########################6 拟时分析(轨迹分析)############################
#准备细胞轨迹分析需要的文件
monocle.matrix=as.matrix(pbmc@assays$RNA@counts)
# monocle.sample=pbmc@meta.data
# monocle.geneAnn=data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
# monocle.clusterAnn=clusterAnn
# monocle.markers=sig.markers

#将Seurat结果转换为monocle需要的细胞矩阵，细胞注释表格和基因注释表格
data <- as(as.matrix(monocle.matrix), 'sparseMatrix')
pd<-new("AnnotatedDataFrame", data = pbmc@meta.data)
fData <- data.frame(gene_short_name = row.names(data),row.names = row.names(data))
fd<-new("AnnotatedDataFrame", data = fData)
cds <- newCellDataSet(data, phenoData = pd, featureData = fd,
                      expressionFamily = negbinomial.size())

# names(pData(cds))[names(pData(cds))=="seurat_clusters"]="Cluster"
# pData(cds)[,"Cluster"]=paste0("cluster",pData(cds)[,"Cluster"])
# 
# #添加细胞聚类数据
# clusterAnn=as.character(monocle.clusterAnn[,2])
# names(clusterAnn)=paste0("cluster",monocle.clusterAnn[,1])
# pData(cds)$cell_type2 <- plyr::revalue(as.character(pData(cds)$Cluster),clusterAnn)

#细胞轨迹分析流程
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds,min_expr = 0.1)

disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table,mean_expression >= 0.1)
cds_myo <- setOrderingFilter(cds,unsup_clustering_genes$gene_id)
plot_ordering_genes(cds_myo)

cds_myo <- reduceDimension(cds_myo,max_components = 2,reduction_method = "DDRTree")
cds_myo <- orderCells(cds_myo) # 包的错 重新安装 解包然后安装知乎的方法删除 
# cds <- setOrderingFilter(cds, as.vector(sig.markers$gene))
# #plot_ordering_genes(cds)
# cds <- reduceDimension(cds, max_components = 2, reduction_method = 'DDRTree')
# cds <- orderCells(cds)

# ######################
# 修改好的包运行了orderCells就可以继续运行下面的结果啦
# ######################
# #保存树枝的细胞轨迹图
# # pdf(file="05.trajectory.State.pdf",width=6.5,height=6)
# plot_cell_trajectory(cds,color_by = "State")
# # dev.off()
# #保存时间的细胞轨迹图
# pdf(file="05.trajectory.Pseudotime.pdf",width=6.5,height=6)
# plot_cell_trajectory(cds,color_by = "Pseudotime")
# dev.off()
# #保存细胞名称的细胞轨迹图
# pdf(file="05.trajectory.cellType.pdf",width=6.5,height=6)
# plot_cell_trajectory(cds,color_by = "cell_type2")
# dev.off()
# #保存聚类的细胞轨迹图
# pdf(file="05.trajectory.cluster.pdf",width=6.5,height=6)
# plot_cell_trajectory(cds, color_by = "Cluster")
# dev.off()
# 
# #细胞轨迹差异分析
# groups=subset(pData(cds),select='State')
# pbmc=AddMetaData(object=pbmc, metadata=groups, col.name="group")
# geneList=list()
# for(i in levels(factor(groups$State))){
#   pbmc.markers=FindMarkers(pbmc, ident.1 = i, group.by = 'group')
#   sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
#   sig.markers=cbind(Gene=row.names(sig.markers), sig.markers)
#   write.table(sig.markers,file=paste0("05.monocleDiff.", i, ".txt"),sep="\t",row.names=F,quote=F)
#   geneList[[i]]=row.names(sig.markers)
# }
# #保存交集基因
# unionGenes=Reduce(union,geneList)
# write.table(file="05.monocleDiff.union.txt",unionGenes,sep="\t",quote=F,col.names=F,row.names=F)





