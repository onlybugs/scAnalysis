######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("Seurat")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("GSVA")
#BiocManager::install("GSEABase")
#BiocManager::install("limma")
# BiocManager::install("SingleR")
# BiocManager::install("celldex")
# BiocManager::install("monocle")



###################################01.数据前期处理和矫正###################################
#引用包
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(celldex)
library(SingleR)
library(monocle)

logFCfilter=1               #logFC的过滤条件
adjPvalFilter=0.05          #矫正后的pvalue的过滤条件
inputFile="treat.txt"       #输入文件



#读取文件，并对输入文件进行整理
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#将矩阵转换为Seurat对象，并对数据进行过滤
pbmc=CreateSeuratObject(counts = data,project = "seurat", min.cells=3, min.features=50, names.delim = "_")
#使用PercentageFeatureSet函数计算线粒体基因的百分比
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
#绘制基因特征的小提琴图
pdf(file="01.featureViolin.pdf", width=10, height=6)
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
pbmc=subset(x = pbmc, subset = nFeature_RNA > 50 & percent.mt < 5)    #对数据进行过滤

#测序深度的相关性图
pdf(file="01.featureCor.pdf",width=10,height=6)
plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5)
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size=1.5)
CombinePlots(plots = list(plot1, plot2))
dev.off()

#对数据进行标准化
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#提取细胞间变异系数较大的基因
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 1500)
#输出特征方差图
top10 <- head(x = VariableFeatures(object = pbmc), 10)
pdf(file="01.featureVar.pdf",width=10,height=6)
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()





###################################02.PCA主成分分析###################################
##PCA分析
pbmc=ScaleData(pbmc)          #PCA降维之前的标准预处理步骤
pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc))     #PCA分析

#绘制每个PCA成分的特征基因
pdf(file="02.pcaGene.pdf",width=10,height=8)
VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca",nfeatures = 20)
dev.off()

#绘制主成分分析图形
pdf(file="02.PCA.pdf",width=6.5,height=6)
DimPlot(object = pbmc, reduction = "pca")
dev.off()

#主成分分析热图
pdf(file="02.pcaHeatmap.pdf",width=10,height=8)
DimHeatmap(object = pbmc, dims = 1:4, cells = 500, balanced = TRUE,nfeatures = 30,ncol=2)
dev.off()

#得到每个PC的p值分布
pbmc <- JackStraw(object = pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:15)
pdf(file="02.pcaJackStraw.pdf",width=8,height=6)
JackStrawPlot(object = pbmc, dims = 1:15)
dev.off()





###################################03.TSNE聚类分析和marker基因###################################
##TSNE聚类分析
pcSelect=14
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)       #计算邻接距离
pbmc <- FindClusters(object = pbmc, resolution = 0.5)         #对细胞分组,对细胞标准模块化
pbmc <- RunTSNE(object = pbmc, dims = 1:pcSelect)             #TSNE聚类
pdf(file="03.TSNE.pdf",width=6.5,height=6)
TSNEPlot(object = pbmc, pt.size = 2, label = TRUE)    #TSNE可视化
dev.off()
write.table(pbmc$seurat_clusters,file="03.tsneCluster.txt",quote=F,sep="\t",col.names=F)

##查找每个聚类的差异基因
pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter)
sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.markers,file="03.clusterMarkers.txt",sep="\t",row.names=F,quote=F)

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#绘制marker在各个cluster的热图
pdf(file="03.tsneHeatmap.pdf",width=12,height=9)
DoHeatmap(object = pbmc, features = top10$gene) + NoLegend()
dev.off()

#绘制marker的小提琴图
pdf(file="03.markerViolin.pdf",width=10,height=6)
VlnPlot(object = pbmc, features = row.names(sig.markers)[1:2])
dev.off()

#需要展示的基因，可以修改
showGenes=c("OLFM4","PDIA2","CPS1","LGALS2","TMED3","CLDN2","MSMB","AQP5","SPINK4","PIGR") 

#绘制marker在各个cluster的散点图
pdf(file="03.markerScatter.pdf",width=10,height=6)
FeaturePlot(object = pbmc, features = showGenes, cols = c("green", "red"))
dev.off()

#绘制marker在各个cluster的气泡图
pdf(file="03.markerBubble.pdf",width=12,height=6)
cluster10Marker=showGenes
DotPlot(object = pbmc, features = cluster10Marker)
dev.off()




###################################04.SingleR R包注释细胞类型###################################
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
write.table(clusterAnn,file="04.clusterAnn.txt",quote=F,sep="\t", row.names=F)
singler2=SingleR(test=counts, ref =ref, 
                labels=ref$label.main)
cellAnn=as.data.frame(singler2)
cellAnn=cbind(id=row.names(cellAnn), cellAnn)
cellAnn=cellAnn[,c("id", "labels")]
write.table(cellAnn, file="04.cellAnn.txt", quote=F, sep="\t", row.names=F)

#cluster注释后的可视化
newLabels=singler$labels
names(newLabels)=levels(pbmc)
pbmc=RenameIdents(pbmc, newLabels)
pdf(file="04.TSNE.pdf",width=6.5,height=6)
TSNEPlot(object = pbmc, pt.size = 2, label = TRUE)    #TSNE可视化
dev.off()

##cluster注释后的差异分析
pbmc.markers=FindAllMarkers(object = pbmc,
                            only.pos = FALSE,
                            min.pct = 0.25,
                            logfc.threshold = logFCfilter)
sig.cellMarkers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.cellMarkers,file="04.cellMarkers.txt",sep="\t",row.names=F,quote=F)



###################################05.monocle R包细胞轨迹分析###################################
#准备细胞轨迹分析需要的文件
monocle.matrix=as.matrix(pbmc@assays$RNA@data)
monocle.sample=pbmc@meta.data
monocle.geneAnn=data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
monocle.clusterAnn=clusterAnn
monocle.markers=sig.markers

#将Seurat结果转换为monocle需要的细胞矩阵，细胞注释表格和基因注释表格
data <- as(as.matrix(monocle.matrix), 'sparseMatrix')
pd<-new("AnnotatedDataFrame", data = monocle.sample)
fd<-new("AnnotatedDataFrame", data = monocle.geneAnn)
cds <- newCellDataSet(data, phenoData = pd, featureData = fd)
names(pData(cds))[names(pData(cds))=="seurat_clusters"]="Cluster"
pData(cds)[,"Cluster"]=paste0("cluster",pData(cds)[,"Cluster"])

#添加细胞聚类数据
clusterAnn=as.character(monocle.clusterAnn[,2])
names(clusterAnn)=paste0("cluster",monocle.clusterAnn[,1])
pData(cds)$cell_type2 <- plyr::revalue(as.character(pData(cds)$Cluster),clusterAnn)

#细胞轨迹分析流程
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- setOrderingFilter(cds, as.vector(sig.markers$gene))
#plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2, reduction_method = 'DDRTree')
cds <- orderCells(cds)
#保存树枝的细胞轨迹图
pdf(file="05.trajectory.State.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "State")
dev.off()
#保存时间的细胞轨迹图
pdf(file="05.trajectory.Pseudotime.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "Pseudotime")
dev.off()
#保存细胞名称的细胞轨迹图
pdf(file="05.trajectory.cellType.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "cell_type2")
dev.off()
#保存聚类的细胞轨迹图
pdf(file="05.trajectory.cluster.pdf",width=6.5,height=6)
plot_cell_trajectory(cds, color_by = "Cluster")
dev.off()

#细胞轨迹差异分析
groups=subset(pData(cds),select='State')
pbmc=AddMetaData(object=pbmc, metadata=groups, col.name="group")
geneList=list()
for(i in levels(factor(groups$State))){
	pbmc.markers=FindMarkers(pbmc, ident.1 = i, group.by = 'group')
	sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
	sig.markers=cbind(Gene=row.names(sig.markers), sig.markers)
	write.table(sig.markers,file=paste0("05.monocleDiff.", i, ".txt"),sep="\t",row.names=F,quote=F)
	geneList[[i]]=row.names(sig.markers)
}
#保存交集基因
unionGenes=Reduce(union,geneList)
write.table(file="05.monocleDiff.union.txt",unionGenes,sep="\t",quote=F,col.names=F,row.names=F)


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
