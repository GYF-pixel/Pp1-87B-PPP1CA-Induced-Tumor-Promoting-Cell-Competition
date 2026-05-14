#Step 1 整合不同Genotype的数据
library(devtools)
library(DropletUtils)
library(Matrix)
library(irlba)
library(Seurat)
library(scran)
library(scDblFinder)
library(scater)
library(DoubletFinder)
library(ggplot2)
library(ggsci)
library(cowplot)
library(tidyverse)
library(ggunchull)
library(SCENIC)
library(scales)
set.seed(1234)
getwd()
setwd("I:\\PPP1CA&PPP1CC\\03scRNA-seq_Mice\\01Analysis\\4th\\Human")

#CA Group
pbmc_CA <- readRDS("CCA_CA.Rds")
pbmc <- pbmc_CA 

##检查default assay
DefaultAssay(pbmc)

##查看有多少个cluster
levels(pbmc)
DefaultAssay(pbmc) <- "RNA"

#提取GFP 细胞
table_GFP <- as.data.frame(pbmc@assays$RNA@counts["GFP",])
table(table_GFP[,1] != 0)
a_GFP <- subset(table_GFP, table_GFP[,1]!=0)

scRNAsub_GFP <- subset(pbmc, cells=row.names(a_GFP))
table(scRNAsub_GFP@meta.data$group)

#提取mCherry 细胞
table_mCherry <- as.data.frame(pbmc@assays$RNA@counts["mCherry",])
table(table_mCherry[,1] != 0)
a_mCherry <- subset(table_mCherry, table_mCherry[,1]!=0)

scRNAsub_mCherry <- subset(pbmc, cells=row.names(a_mCherry))
table(scRNAsub_mCherry@meta.data$group)

DefaultAssay(scRNAsub_mCherry)
DefaultAssay(scRNAsub_GFP)


CA_scRNAsub_GFP_mCherry.anchors <- FindIntegrationAnchors(object.list = list(scRNAsub_GFP, scRNAsub_mCherry), anchor.features = 2000, dims = 1:50)
CA_scRNAsub_GFP_mCherry.combined <- IntegrateData(anchorset = CA_scRNAsub_GFP_mCherry.anchors, dims = 1:50)
pbmc <-  CA_scRNAsub_GFP_mCherry.combined

#pbmc <- NormalizeData(pbmc, verbose = FALSE, normalization.method = "LogNormalize", scale.factor = 1e4)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2500)
pbmc <- ScaleData(pbmc, vars.to.regress = c("nCount_RNA"), verbose = TRUE)
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc), npcs = 40, nfeature.print = 10, ndims.print = 1:5, verbose = T)
pc.num=1:40
pbmc <- RunUMAP(pbmc, dims=pc.num)
pbmc <- FindNeighbors(pbmc, dims = pc.num)
pbmc = FindClusters(pbmc,resolution = 1.0)

saveRDS(pbmc, file="CCA_CA_GFP_mCherry.Rds") 


################
pbmc <- readRDS("CCA_CA_GFP_mCherry.Rds")

##检查default assay
DefaultAssay(pbmc)

##查看有多少个cluster
levels(pbmc)

##设置对应cluster的配色数目
col5 <- colorRampPalette((pal_npg(palette = c("nrc"))(7)))(11)
show_col(col5)

##正常的UMAP降维数据结果
UMPplot_label <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("Integrated scRNA-seq datasets")
ggsave("P1_UMPplot_label.pdf", plot = plot_grid(UMPplot_label), width = 6, height = 5)
UMPplot_unlabel <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("Integrated scRNA-seq datasets")
ggsave("P1_UMPplot_unlabel.pdf", plot = plot_grid(UMPplot_unlabel), width = 6, height = 5)


##相同Genotype的分组结果
table(pbmc@meta.data[["group"]])
##按group区别
UMPplot_label_group <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "group", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("Integrated scRNA-seq datasets")
ggsave("P2_UMPplot_label_group.pdf", plot = plot_grid(UMPplot_label_group), width = 6, height = 5)
UMPplot_unlabel_group <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "group", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("Integrated scRNA-seq datasets")
ggsave("P2_UMPplot_unlabel_group.pdf", plot = plot_grid(UMPplot_unlabel_group), width = 6, height = 5)

##相同Genotype的分组结果
table(pbmc@meta.data[["sec_group"]])
##按sec_group区别
UMPplot_label_sec_group <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "sec_group", pt.size = 1.0, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("Integrated scRNA-seq datasets")
ggsave("P2_UMPplot_label_sec_group.pdf", plot = plot_grid(UMPplot_label_sec_group), width = 10, height = 8)
UMPplot_unlabel_sec_group <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "sec_group", pt.size = 1.0, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("Integrated scRNA-seq datasets")
ggsave("P2_UMPplot_unlabel_sec_group.pdf", plot = plot_grid(UMPplot_unlabel_sec_group), width = 10, height = 8)

##各cell cluster的细胞占比
table(pbmc@meta.data$group, pbmc@meta.data$integrated_snn_res.1)
table(pbmc@meta.data$sec_group, pbmc@meta.data$integrated_snn_res.1)

##生成文件画出堆积柱状图
Cells_in_different_Clusters_at_different_group <- table(pbmc@meta.data$sec_group, pbmc@meta.data$integrated_snn_res.1)
write.csv(Cells_in_different_Clusters_at_different_group,file = "P3_Cells_in_different_Clusters_at_different_group.csv",row.names = T)
##计算单个cluster中两种cell type的占比，设置Unique Cluster筛选条件


#Marker Gene Expression
a <- FeaturePlot(pbmc,features = c("GFP"),cols = c('lightgray','red'), pt.size = 1.0)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("GFP")
ggsave("P6_pbmcUMAP1_GFPplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("mCherry"),cols = c('lightgray','red'), pt.size = 1.0)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("mCherry")
ggsave("P6_pbmcUMAP1_mCherryplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("PPP1CA"),cols = c('lightgray','red'), pt.size = 1.0)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("PPP1CA")
ggsave("P6_pbmcUMAP1_PPP1CAplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("PPP1CC"),cols = c('lightgray','red'), pt.size = 1.0)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("PPP1CC")
ggsave("P6_pbmcUMAP1_PPP1CCplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("KRAS"),cols = c('lightgray','red'), pt.size = 1.0)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("KRAS")
ggsave("P6_pbmcUMAP1_KRASplot2.pdf", plot = plot_grid(a), width = 6, height = 5)



                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     



###subset cluster 3
#Select Unique Cluster 
pbmc@meta.data$UniqueCluster = NA

Cells.sub <- subset(pbmc@meta.data, integrated_snn_res.1 %in% c("3"))
table(pbmc@meta.data$integrated_snn_res.1)
table(Cells.sub$integrated_snn_res.1)


scRNAsub <- subset(pbmc, cells=row.names(Cells.sub))

##重新聚类分析
DefaultAssay(scRNAsub) <- "RNA"

scRNAsub <- FindVariableFeatures(scRNAsub, selection.method = "vst", nfeatures = 2500)
scRNAsub <- ScaleData(scRNAsub, vars.to.regress = c("nCount_RNA"), verbose = TRUE)
scRNAsub <- RunPCA(scRNAsub, features = VariableFeatures(scRNAsub), npcs = 60, nfeature.print = 10, ndims.print = 1:5, verbose = T)
pc.num=1:60
scRNAsub <- RunUMAP(scRNAsub, dims=pc.num)
scRNAsub <- FindNeighbors(scRNAsub, dims = pc.num)
scRNAsub = FindClusters(scRNAsub,resolution = 1.2)

##查看有多少个cluster
levels(scRNAsub)

##设置对应cluster的配色数目
col5 <- colorRampPalette((pal_npg(palette = c("nrc"))(7)))(15)
show_col(col5)

##正常的UMAP降维数据结果
UMPplot_label <- DimPlot(scRNAsub, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("WT & RasRcd5 Epithelium")
ggsave("P7_UMPplot_label.pdf", plot = plot_grid(UMPplot_label), width = 6, height = 5)
UMPplot_unlabel <- DimPlot(scRNAsub, reduction = "umap", label = FALSE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("WT & RasRcd5 Epithelium")
ggsave("P7_UMPplot_unlabel.pdf", plot = plot_grid(UMPplot_unlabel), width = 6, height = 5)

#Marker Gene Expression
a <- FeaturePlot(scRNAsub,features = c("GFP"),cols = c('lightgray','red'), pt.size = 1.0)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("GFP")
ggsave("P8_scRNAsubUMAP1_GFPplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(scRNAsub,features = c("mCherry"),cols = c('lightgray','red'), pt.size = 1.0)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("mCherry")
ggsave("P8_scRNAsubUMAP1_mCherryplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(scRNAsub,features = c("PPP1CA"),cols = c('lightgray','red'), pt.size = 1.0)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("PPP1CA")
ggsave("P8_scRNAsubUMAP1_PPP1CAplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(scRNAsub,features = c("PPP1CC"),cols = c('lightgray','red'), pt.size = 1.0)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("PPP1CC")
ggsave("P8_scRNAsubUMAP1_PPP1CCplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(scRNAsub,features = c("KRAS"),cols = c('lightgray','red'), pt.size = 1.0)+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("KRAS")
ggsave("P8_scRNAsubUMAP1_KRASplot2.pdf", plot = plot_grid(a), width = 6, height = 5)













#####AddModuleScore
library(RColorBrewer)

gene_sets <- read.table("sgCA_vs_sgNT_Down_genes.txt")
gene_sets2 <- list(gene_sets[,1])
gene_sets2

sce_T <- AddModuleScore(pbmc, features = gene_sets2, nbin = 24, ctrl = 100,name = "sgCA_vs_sgNT_Down_genes")
FeaturePlot(sce_T,'sgCA_vs_sgNT_Down_genes1',cols=rev(brewer.pal(10, name = "RdBu")))
ggsave("P7_CA_sgCA_vs_sgNT_Down_genes.pdf")


gene_sets <- read.table("sgCA_vs_sgNT_Up_genes.txt")
gene_sets2 <- list(gene_sets[,1])
gene_sets2

sce_T <- AddModuleScore(pbmc, features = gene_sets2, nbin = 24, ctrl = 100,name = "sgCA_vs_sgNT_Up_genes")
FeaturePlot(sce_T,'sgCA_vs_sgNT_Up_genes1',cols=rev(brewer.pal(10, name = "RdBu")))
ggsave("P7_CA_sgCA_vs_sgNT_Up_genes.pdf")





