#Step 1 整合不同Genotype的数据
library(devtools)
library(DropletUtils)
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
library(AUCell)
set.seed(1234)
getwd()
setwd("I:\\PPP1CA&PPP1CC\\03scRNA-seq_Mice\\01Analysis\\4th_Try")
setwd("H:\\PPP1CA&PPP1CC\\03scRNA-seq_Mice\\01Analysis\\4th_Try")


pbmc_WT <- readRDS("CCA_NT_over1500.Rds")

#WT Group
pbmc <- pbmc_WT

##检查default assay
DefaultAssay(pbmc)

##查看有多少个cluster
levels(pbmc)

##设置对应cluster的配色数目
col5 <- colorRampPalette((pal_npg(palette = c("nrc"))(7)))(15)
show_col(col5)

##正常的UMAP降维数据结果
UMPplot_label <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("WT scRNA-seq datasets")
ggsave("P1_UMPplot_label.pdf", plot = plot_grid(UMPplot_label), width = 10, height = 8)
UMPplot_unlabel <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("WT scRNA-seq datasets")
ggsave("P1_UMPplot_unlabel.pdf", plot = plot_grid(UMPplot_unlabel), width = 10, height = 8)

##相同Genotype的分组结果
table(pbmc@meta.data[["group"]])
##按group区别
UMPplot_label_group <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "group", pt.size = 1.0, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("WT scRNA-seq datasets")
ggsave("P2_UMPplot_label_group.pdf", plot = plot_grid(UMPplot_label_group), width = 10, height = 8)
UMPplot_unlabel_group <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "group", pt.size = 1.0, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("WT scRNA-seq datasets")
ggsave("P2_UMPplot_unlabel_group.pdf", plot = plot_grid(UMPplot_unlabel_group), width = 10, height = 8)

##各cell cluster的细胞占比
table(pbmc@meta.data$group, pbmc@meta.data$integrated_snn_res.0.7)
##生成文件画出堆积柱状图
Cells_in_different_Clusters_at_different_group <- table(pbmc@meta.data$group, pbmc@meta.data$integrated_snn_res.0.7)
write.csv(Cells_in_different_Clusters_at_different_group,file = "P3_Cells_in_different_Clusters_at_different_group.csv",row.names = T)

##各cell cluster--Find markers for every cluster compared to all remaining cells
##这里需要将默认格式输出调整为“RNA”
DefaultAssay(pbmc) <- "RNA"
pbmc <- ScaleData(pbmc, vars.to.regress = c("nCount_RNA"), verbose = TRUE)
pbmc <- ScaleData(pbmc, vars.to.regress = c("S.Score", "G2M.Score","percent.mt", "nCount_RNA"), verbose = TRUE)

## report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc,
                               only.pos = FALSE,
                               min.pct = 0.4,
                               logfc.threshold = 0.25)

Allmarkers_pbmc = pbmc.markers %>% select(gene, everything()) %>% subset(p_val<0.05)
write.csv(Allmarkers_pbmc, "P4_Allmarkers__pbmc_wilcox.csv", row.names = T)

## 使用Top 5基因画热图
top5pbmc.markers <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
##做Heatmap的热图
Heatmap_markers_pbmc <- DoHeatmap(pbmc,features = top5pbmc.markers$gene,
          group.colors = col5) + theme(text = element_text (size = 18)) + 
		  scale_fill_gradient2(low = '#0099CC',mid = 'white',high = '#CC0033',name = 'Z-score')
ggplot2::ggsave(filename = 'P4_Heatmap_markers_pbmc.pdf',width =16,height = 19)

##做jjVolcano的火山图
library(scRNAtoolVis)
jjVolcano_scRNAsub.markers <- jjVolcano(diffData = pbmc.markers,
          log2FC.cutoff = 0.25, 
          size  = 3.5, #设置点的大小
          fontface = 'italic', #设置字体形式
          aesCol = c('#00468B','#ED0000'), #设置点的颜色
          tile.col = col5, #设置cluster的颜色
          #col.type = "adjustP", #设置矫正方式
          topGeneN = 5 #设置展示topN的基因
         )
jjVolcano_scRNAsub.markers
ggsave(filename = 'P4_jjVolcano_scRNAsub.markers.pdf', plot = plot_grid(jjVolcano_scRNAsub.markers), width = 15, height = 6)

#Marker Gene Expression
a <- FeaturePlot(pbmc,features = c("GFP"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("GFP")
ggsave("P6_pbmcUMAP1_GFPplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("mCherry"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("mCherry")
ggsave("P6_pbmcUMAP1_mCherryplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("PPP1CA"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("PPP1CA")
ggsave("P6_pbmcUMAP1_PPP1CAplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("PPP1CC"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("PPP1CC")
ggsave("P6_pbmcUMAP1_PPP1CCplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("KRAS"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("KRAS")
ggsave("P6_pbmcUMAP1_KRASplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("Cas9"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Cas9")
ggsave("P6_pbmcUMAP1_Cas9plot2.pdf", plot = plot_grid(a), width = 6, height = 5)

#提取GFP 细胞
table_GFP <- as.data.frame(pbmc@assays$RNA@counts["GFP",])
table(table_GFP[,1] != 0)
a_GFP <- subset(table_GFP, table_GFP[,1]!=0)

scRNAsub_GFP <- subset(pbmc, cells=row.names(a_GFP))
table(scRNAsub_GFP@meta.data$group)

GFP_seurat_clusters <- table(scRNAsub_GFP@meta.data$sec_group, scRNAsub_GFP@meta.data$seurat_clusters)
GFP_seurat_clusters 
write.csv(GFP_seurat_clusters, file = "P4_GFP_seurat_clusters_WT.csv",row.names = T)
saveRDS(scRNAsub_GFP, file="CCA_NT_scRNAsub_GFP.Rds") 

#提取mCherry 细胞
table_mCherry <- as.data.frame(pbmc@assays$RNA@counts["mCherry",])
table(table_mCherry[,1] != 0)
a_mCherry <- subset(table_mCherry, table_mCherry[,1]!=0)

scRNAsub_mCherry <- subset(pbmc, cells=row.names(a_mCherry))
table(scRNAsub_mCherry@meta.data$group)

mCherry_seurat_clusters <- table(scRNAsub_mCherry@meta.data$sec_group, scRNAsub_mCherry@meta.data$seurat_clusters)
mCherry_seurat_clusters 
write.csv(mCherry_seurat_clusters, file = "P4_mCherry_seurat_clusters_WT.csv",row.names = T)
saveRDS(scRNAsub_mCherry, file="CCA_NT_scRNAsub_mCherry.Rds") 



#CA Group
pbmc_CA <- readRDS("CCA_CA.Rds")
pbmc <- pbmc_CA 

##检查default assay
DefaultAssay(pbmc)

##查看有多少个cluster
levels(pbmc)

##设置对应cluster的配色数目
col5 <- colorRampPalette((pal_npg(palette = c("nrc"))(7)))(14)
show_col(col5)

##正常的UMAP降维数据结果
UMPplot_label <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("sgPPP1CA scRNA-seq datasets")
ggsave("P1_UMPplot_label.pdf", plot = plot_grid(UMPplot_label), width = 10, height = 8)
UMPplot_unlabel <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("sgPPP1CA scRNA-seq datasets")
ggsave("P1_UMPplot_unlabel.pdf", plot = plot_grid(UMPplot_unlabel), width = 10, height = 8)

##相同Genotype的分组结果
table(pbmc@meta.data[["group"]])
##按group区别
UMPplot_label_group <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "group", pt.size = 1.0, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("sgPPP1CA scRNA-seq datasets")
ggsave("P2_UMPplot_label_group.pdf", plot = plot_grid(UMPplot_label_group), width = 10, height = 8)
UMPplot_unlabel_group <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "group", pt.size = 1.0, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("sgPPP1CA scRNA-seq datasets")
ggsave("P2_UMPplot_unlabel_group.pdf", plot = plot_grid(UMPplot_unlabel_group), width = 10, height = 8)

##按phase区别
UMPplot_unlabel_group <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "Phase", pt.size = 1.0, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("sgPPP1CA scRNA-seq datasets")
ggsave("P2_UMPplot_unlabel_Phase.pdf", plot = plot_grid(UMPplot_unlabel_group), width = 10, height = 8)

##各cell cluster的细胞占比
table(pbmc@meta.data$group, pbmc@meta.data$RNA_snn_res.0.9)
##生成文件画出堆积柱状图
Cells_in_different_Clusters_at_different_group <- table(pbmc@meta.data$group, pbmc@meta.data$RNA_snn_res.0.9)
write.csv(Cells_in_different_Clusters_at_different_group,file = "P3_Cells_in_different_Clusters_at_different_group.csv",row.names = T)

## report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc,
                               only.pos = FALSE,
                               min.pct = 0.4,
                               logfc.threshold = 0.25)

Allmarkers_pbmc = pbmc.markers %>% select(gene, everything()) %>% subset(p_val<0.05)
write.csv(Allmarkers_pbmc, "P4_Allmarkers__pbmc_wilcox.csv", row.names = T)

## 使用Top 5基因画热图
top5pbmc.markers <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
##做Heatmap的热图
Heatmap_markers_pbmc <- DoHeatmap(pbmc,features = top5pbmc.markers$gene,
          group.colors = col5) + theme(text = element_text (size = 18)) + 
		  scale_fill_gradient2(low = '#0099CC',mid = 'white',high = '#CC0033',name = 'Z-score')
ggplot2::ggsave(filename = 'P4_Heatmap_markers_pbmc.pdf',width =16,height = 19)

##做jjVolcano的火山图
library(scRNAtoolVis)
jjVolcano_scRNAsub.markers <- jjVolcano(diffData = pbmc.markers,
          log2FC.cutoff = 0.25, 
          size  = 3.5, #设置点的大小
          fontface = 'italic', #设置字体形式
          aesCol = c('#00468B','#ED0000'), #设置点的颜色
          tile.col = col5, #设置cluster的颜色
          #col.type = "adjustP", #设置矫正方式
          topGeneN = 5 #设置展示topN的基因
         )
jjVolcano_scRNAsub.markers
ggsave(filename = 'P4_jjVolcano_scRNAsub.markers.pdf', plot = plot_grid(jjVolcano_scRNAsub.markers), width = 15, height = 6)

#Marker Gene Expression
a <- FeaturePlot(pbmc,features = c("GFP"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("GFP")
ggsave("P6_pbmcUMAP1_GFPplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("mCherry"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("mCherry")
ggsave("P6_pbmcUMAP1_mCherryplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("PPP1CA"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("PPP1CA")
ggsave("P6_pbmcUMAP1_PPP1CAplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("PPP1CC"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("PPP1CC")
ggsave("P6_pbmcUMAP1_PPP1CCplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("KRAS"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("KRAS")
ggsave("P6_pbmcUMAP1_KRASplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("Cas9"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Cas9")
ggsave("P6_pbmcUMAP1_Cas9plot2.pdf", plot = plot_grid(a), width = 6, height = 5)











#CA Group
pbmc_CA <- readRDS("CCA_CA_over1500genes.Rds")
pbmc <- pbmc_CA 

##检查default assay
DefaultAssay(pbmc)

##查看有多少个cluster
levels(pbmc)

##设置对应cluster的配色数目
col5 <- colorRampPalette((pal_npg(palette = c("nrc"))(7)))(13)
show_col(col5)

##正常的UMAP降维数据结果
UMPplot_label <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("sgPPP1CA scRNA-seq datasets")
ggsave("P1_UMPplot_label.pdf", plot = plot_grid(UMPplot_label), width = 10, height = 8)
UMPplot_unlabel <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("sgPPP1CA scRNA-seq datasets")
ggsave("P1_UMPplot_unlabel.pdf", plot = plot_grid(UMPplot_unlabel), width = 10, height = 8)

##相同Genotype的分组结果
table(pbmc@meta.data[["group"]])
##按group区别
UMPplot_label_group <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "group", pt.size = 1.0, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("sgPPP1CA scRNA-seq datasets")
ggsave("P2_UMPplot_label_group.pdf", plot = plot_grid(UMPplot_label_group), width = 10, height = 8)
UMPplot_unlabel_group <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "group", pt.size = 1.0, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("sgPPP1CA scRNA-seq datasets")
ggsave("P2_UMPplot_unlabel_group.pdf", plot = plot_grid(UMPplot_unlabel_group), width = 10, height = 8)

##按phase区别
UMPplot_unlabel_group <- DimPlot(pbmc, reduction = "umap", label = FALSE, group.by = "Phase", pt.size = 1.0, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("sgPPP1CA scRNA-seq datasets")
ggsave("P2_UMPplot_unlabel_Phase.pdf", plot = plot_grid(UMPplot_unlabel_group), width = 10, height = 8)

##各cell cluster的细胞占比
table(pbmc@meta.data$group, pbmc@meta.data$integrated_snn_res.0.7)
##生成文件画出堆积柱状图
Cells_in_different_Clusters_at_different_group <- table(pbmc@meta.data$group, pbmc@meta.data$integrated_snn_res.0.7)
write.csv(Cells_in_different_Clusters_at_different_group,file = "P3_Cells_in_different_Clusters_at_different_group.csv",row.names = T)

##各cell cluster--Find markers for every cluster compared to all remaining cells
##这里需要将默认格式输出调整为“RNA”
DefaultAssay(pbmc) <- "RNA"
pbmc <- ScaleData(pbmc, vars.to.regress = c("nCount_RNA"), verbose = TRUE)
pbmc <- ScaleData(pbmc, vars.to.regress = c("S.Score", "G2M.Score","percent.mt", "nCount_RNA"), verbose = TRUE)

## report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc,
                               only.pos = FALSE,
                               min.pct = 0.4,
                               logfc.threshold = 0.25)

Allmarkers_pbmc = pbmc.markers %>% select(gene, everything()) %>% subset(p_val<0.05)
write.csv(Allmarkers_pbmc, "P4_Allmarkers__pbmc_wilcox.csv", row.names = T)

## 使用Top 5基因画热图
top5pbmc.markers <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
##做Heatmap的热图
Heatmap_markers_pbmc <- DoHeatmap(pbmc,features = top5pbmc.markers$gene,
          group.colors = col5) + theme(text = element_text (size = 18)) + 
		  scale_fill_gradient2(low = '#0099CC',mid = 'white',high = '#CC0033',name = 'Z-score')
ggplot2::ggsave(filename = 'P4_Heatmap_markers_pbmc.pdf',width =16,height = 19)

##做jjVolcano的火山图
library(scRNAtoolVis)
jjVolcano_scRNAsub.markers <- jjVolcano(diffData = pbmc.markers,
          log2FC.cutoff = 0.25, 
          size  = 3.5, #设置点的大小
          fontface = 'italic', #设置字体形式
          aesCol = c('#00468B','#ED0000'), #设置点的颜色
          tile.col = col5, #设置cluster的颜色
          #col.type = "adjustP", #设置矫正方式
          topGeneN = 5 #设置展示topN的基因
         )
jjVolcano_scRNAsub.markers
ggsave(filename = 'P4_jjVolcano_scRNAsub.markers.pdf', plot = plot_grid(jjVolcano_scRNAsub.markers), width = 15, height = 6)

#Marker Gene Expression
a <- FeaturePlot(pbmc,features = c("GFP"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("GFP")
ggsave("P6_pbmcUMAP1_GFPplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("mCherry"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("mCherry")
ggsave("P6_pbmcUMAP1_mCherryplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("PPP1CA"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("PPP1CA")
ggsave("P6_pbmcUMAP1_PPP1CAplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("PPP1CC"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("PPP1CC")
ggsave("P6_pbmcUMAP1_PPP1CCplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("KRAS"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("KRAS")
ggsave("P6_pbmcUMAP1_KRASplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(pbmc,features = c("Cas9"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Cas9")
ggsave("P6_pbmcUMAP1_Cas9plot2.pdf", plot = plot_grid(a), width = 6, height = 5)



#提取GFP 细胞
table_GFP <- as.data.frame(pbmc@assays$RNA@counts["GFP",])
table(table_GFP[,1] != 0)
a_GFP <- subset(table_GFP, table_GFP[,1]!=0)

scRNAsub_GFP <- subset(pbmc, cells=row.names(a_GFP))
table(scRNAsub_GFP@meta.data$group)

GFP_seurat_clusters <- table(scRNAsub_GFP@meta.data$sec_group, scRNAsub_GFP@meta.data$seurat_clusters)
GFP_seurat_clusters 
write.csv(GFP_seurat_clusters, file = "P4_GFP_seurat_clusters_CA.csv",row.names = T)
saveRDS(scRNAsub_GFP, file="CCA_CA_scRNAsub_GFP.Rds") 

#提取mCherry 细胞
table_mCherry <- as.data.frame(pbmc@assays$RNA@counts["mCherry",])
table(table_mCherry[,1] != 0)
a_mCherry <- subset(table_mCherry, table_mCherry[,1]!=0)

scRNAsub_mCherry <- subset(pbmc, cells=row.names(a_mCherry))
table(scRNAsub_mCherry@meta.data$group)

mCherry_seurat_clusters <- table(scRNAsub_mCherry@meta.data$sec_group, scRNAsub_mCherry@meta.data$seurat_clusters)
mCherry_seurat_clusters 
write.csv(mCherry_seurat_clusters, file = "P4_mCherry_seurat_clusters_CA.csv",row.names = T)
saveRDS(scRNAsub_mCherry, file="CCA_CA_scRNAsub_mCherry.Rds") 



############################################################

#只用CA1
meta <- pbmc@meta.data
meta2 <- meta[pbmc@meta.data$group == "CA1",]
View(meta2)
scRNAsub2 <- subset(pbmc, cells=row.names(meta2))
table(meta2$group)

#DefaultAssay(scRNAsub2) <- "integrated"
DefaultAssay(scRNAsub2) <- "RNA"
scRNAsub2 <- NormalizeData(scRNAsub2, verbose = FALSE, normalization.method = "LogNormalize", scale.factor = 1e4)
scRNAsub2 <- FindVariableFeatures(scRNAsub2, selection.method = "vst", nfeatures = 2500)
scRNAsub2 <- ScaleData(scRNAsub2, vars.to.regress = c("S.Score", "G2M.Score","percent.mt", "nCount_RNA"), verbose = TRUE)
scRNAsub2 <- RunPCA(scRNAsub2, npcs = 60, verbose = FALSE)
scRNAsub2 <- RunUMAP(scRNAsub2, reduction = "pca", dims = 1:60)
scRNAsub2 <- FindNeighbors(scRNAsub2, dims = 1:60)
scRNAsub2 = FindClusters(scRNAsub2,resolution = 0.9)


##检查default assay
DefaultAssay(scRNAsub2)
##查看有多少个cluster
levels(scRNAsub2)

##设置对应cluster的配色数目
col5 <- colorRampPalette((pal_npg(palette = c("nrc"))(7)))(13)
show_col(col5)

##正常的UMAP降维数据结果
UMPplot_label <- DimPlot(scRNAsub2, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("sgPPP1CA scRNA-seq datasets")
ggsave("P1_UMPplot_label.pdf", plot = plot_grid(UMPplot_label), width = 10, height = 8)
UMPplot_unlabel <- DimPlot(scRNAsub2, reduction = "umap", label = FALSE, group.by = "seurat_clusters", pt.size = 1.0, label.size = 7, cols = col5, raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("sgPPP1CA scRNA-seq datasets")
ggsave("P1_UMPplot_unlabel.pdf", plot = plot_grid(UMPplot_unlabel), width = 10, height = 8)

##相同Genotype的分组结果
table(scRNAsub2@meta.data[["group"]])
##按group区别
UMPplot_label_group <- DimPlot(scRNAsub2, reduction = "umap", label = TRUE, group.by = "group", pt.size = 1.0, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("sgPPP1CA scRNA-seq datasets")
ggsave("P2_UMPplot_label_group.pdf", plot = plot_grid(UMPplot_label_group), width = 10, height = 8)
UMPplot_unlabel_group <- DimPlot(scRNAsub2, reduction = "umap", label = FALSE, group.by = "group", pt.size = 1.0, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("sgPPP1CA scRNA-seq datasets")
ggsave("P2_UMPplot_unlabel_group.pdf", plot = plot_grid(UMPplot_unlabel_group), width = 10, height = 8)

##按phase区别
UMPplot_unlabel_group <- DimPlot(scRNAsub2, reduction = "umap", label = FALSE, group.by = "Phase", pt.size = 1.0, label.size = 7, cols = c("#4DBBD5","#E64B35","#00A087","#3C5488"), raster=FALSE)+theme(panel.border = element_rect(color = "black",linewidth = 2), text = element_text (size = 18),axis.text = element_text (size = 18))+ggtitle("sgPPP1CA scRNA-seq datasets")
ggsave("P2_UMPplot_unlabel_Phase.pdf", plot = plot_grid(UMPplot_unlabel_group), width = 10, height = 8)

##各cell cluster的细胞占比
table(scRNAsub2@meta.data$group, scRNAsub2@meta.data$integrated_snn_res.0.7)
##生成文件画出堆积柱状图
Cells_in_different_Clusters_at_different_group <- table(scRNAsub2@meta.data$group, scRNAsub2@meta.data$integrated_snn_res.0.7)
write.csv(Cells_in_different_Clusters_at_different_group,file = "P3_Cells_in_different_Clusters_at_different_group.csv",row.names = T)

#Marker Gene Expression
a <- FeaturePlot(scRNAsub2,features = c("GFP"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("GFP")
ggsave("P6_scRNAsub2UMAP1_GFPplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(scRNAsub2,features = c("mCherry"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("mCherry")
ggsave("P6_scRNAsub2UMAP1_mCherryplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(scRNAsub2,features = c("PPP1CA"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("PPP1CA")
ggsave("P6_scRNAsub2UMAP1_PPP1CAplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(scRNAsub2,features = c("PPP1CC"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("PPP1CC")
ggsave("P6_scRNAsub2UMAP1_PPP1CCplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(scRNAsub2,features = c("KRAS"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("KRAS")
ggsave("P6_scRNAsub2UMAP1_KRASplot2.pdf", plot = plot_grid(a), width = 6, height = 5)
a <- FeaturePlot(scRNAsub2,features = c("Cas9"),cols = c('lightgray','red'))+theme(panel.border = element_rect(color = "black",linewidth = 2))+ggtitle("Cas9")
ggsave("P6_scRNAsub2UMAP1_Cas9plot2.pdf", plot = plot_grid(a), width = 6, height = 5)
















saveRDS(pbmc, file="CCA_CA.Rds") 

#####AddModuleScore

library(RColorBrewer)
pbmc <- readRDS(file="CCA_CA.Rds")

gene_sets <- read.table("sgCA_vs_sgNT_Down_genes.txt")
gene_sets2 <- list(gene_sets[,1])
gene_sets2

sce_T <- AddModuleScore(pbmc, features = gene_sets2, nbin = 24, ctrl = 100,name = "sgCA_vs_sgNT_Down_genes")
FeaturePlot(sce_T,'sgCA_vs_sgNT_Down_genes1',cols=rev(brewer.pal(10, name = "RdBu")))
ggsave("P6_CA_sgCA_vs_sgNT_Down_genes.pdf")


gene_sets <- read.table("sgCA_vs_sgNT_Up_genes.txt")
gene_sets2 <- list(gene_sets[,1])
gene_sets2

sce_T <- AddModuleScore(pbmc, features = gene_sets2, nbin = 24, ctrl = 100,name = "sgCA_vs_sgNT_Up_genes")
FeaturePlot(sce_T,'sgCA_vs_sgNT_Up_genes1',cols=rev(brewer.pal(10, name = "RdBu")))
ggsave("P6_CA_sgCA_vs_sgNT_Up_genes.pdf")







