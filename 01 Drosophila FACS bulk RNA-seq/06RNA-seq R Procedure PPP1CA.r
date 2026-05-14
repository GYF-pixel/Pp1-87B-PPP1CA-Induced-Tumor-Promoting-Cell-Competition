library(devtools)
library(biomaRt)
library(curl)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(stringr)
library(edgeR)
library(ggsci)
library(cowplot)
library(tidyverse)
library(ggunchull)
library(SCENIC)
library(scales)
#library(xlsx)
setwd("H:\\PPP1CA&PPP1CC\\01BulkRNAseq\\02Analysis\\00CountFile")                    #设置工作目录



##foldChange=0.5849625    foldChange=log2(FC=1.5)=0.5849625
foldChange=0.5849625
padj=0.05

#Step 1 
#整理RNA-seq的基因表达矩阵
#https://mp.weixin.qq.com/s/qsBOPoJwItO5qAGstKRC_A

#创建一个文件夹，这个folder里只有sample的count文件
#显示目前工作目录有sample文件
#all_files <- list.files(path='K:\\E75\\RNA_seq part\\04Expression\\EyeDisc')
all_files <- list.files(path = )
print(all_files)

#提取sampleID
filesID <- as.data.frame(all_files)
print(filesID)
sampleID <- separate(data = filesID,col = all_files, into = c("id","other"),sep = "_gene")
print(sampleID)

#读取并合并不同sample的基因表达矩阵
#显示所有被读取的文件的路径
for (i in 1:length(all_files))
{
  files <- all_files[i]
  model_path = paste('H:\\PPP1CA&PPP1CC\\01BulkRNAseq\\02Analysis\\00CountFile\\',files, sep='')
  print(model_path)
}

#将所有的sample储存在一个list对象里
q <- list()
for (i in 1:length(all_files))
{
  files <- all_files[i]
  model_path = paste('H:\\PPP1CA&PPP1CC\\01BulkRNAseq\\02Analysis\\00CountFile\\',files, sep='')
  data_table = as.data.frame(read.table(model_path,sep = "\t",col.names = c("gene_id", sampleID$id[i])))
  rownames(data_table) <- data_table[,1]
  q[[i]] <- data_table
}

#将所有的sample合并到一个data.frame的对象里
newdata <- as.data.frame(q[[1]])
for (i in 2:length(q))
{
  newdata <- cbind(newdata, q[[i]][2])
}

#查看FBgn和FBti的分割线
newdata[17865:17880,]
newdata[17869:17880,]
newdata_filter <- newdata[1:17869,]
tail(newdata_filter)
newdata_filter <- newdata_filter[,-1]
write.csv(newdata_filter,file ="newdata_filter.csv")

##Step 2 DEGs identification by edgeR
library(devtools)
library(biomaRt)
library(curl)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(stringr)
library(edgeR)
library(ggsci)
library(cowplot)
library(tidyverse)
library(ggunchull)
library(SCENIC)
library(scales)
setwd("H:\\PPP1CA&PPP1CC\\01BulkRNAseq\\02Analysis")                    #设置工作目录
setwd("I:\\PPP1CA&PPP1CC\\01BulkRNAseq\\02Analysis")   

##foldChange=0.5849625    foldChange=log2(FC=1.5)=0.5849625
foldChange=0.5849625
padj=0.05

##手动删除all_gene_count.tsv没有数据对应的行 
rt=read.csv("newdata_filter.csv")              #改成自己的文件名
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
data = exp
data.filter = data[rowMeans(data) > 1 & rowSums(data) > 15, ]   #数据过滤

#####检测Batch Effect 批次效应
#bulk RNA-Seq（6）样本相关性、聚类、PCA分析：https://mp.weixin.qq.com/s/kOpEVIxThhmTEogM89OGLg
#RNAseq常规主成分分析（PCA）不好看，可以这样做: https://mp.weixin.qq.com/s/ja90VBIycMAr9L7aoCtNyA

#1.相关性分析
#计算距离
sample_cor <- cor(data.filter)
sample_cor1 <- round(sample_cor, digits = 2)
#画图
library(pheatmap)
cor_heatmap <- pheatmap(sample_cor1, display_numbers = T,fontsize = 10, angle_col = 45)
cor_heatmap
ggsave("00_BatchEffect_Correlation_heatmap.pdf", plot = cor_heatmap, width = 11, height = 8) 


#2. PCA 分析
gene <- t(data.filter)
library(FactoMineR)
gene.pca <- PCA(gene, ncp = 2, scale.unit = TRUE, graph = FALSE)
plot(gene.pca)  #PCA 简图

#提取样本在 PCA 前两轴中的坐标
pca_sample <- data.frame(gene.pca$ind$coord[ ,1:2])
head(pca_sample)
pca_sample

#提取 PCA 前两轴的贡献度
pca_eig1 <- round(gene.pca$eig[1,2], 2)
pca_eig2 <- round(gene.pca$eig[2,2],2 )
pca_eig1
pca_eig2

#增添分组信息
colnames(data.filter)
group_info = c(rep("Ras_GFP_neg",3), rep("Ras_GFP_pos",3), rep("RasPPP1CA_GFP_neg",3), 
rep("RasPPP1CA_GFP_pos",4), rep("WT_GFP_neg",3), rep("WT_GFP_pos",3))

group_info2 <- as.data.frame(group_info)
rownames(group_info2) <- colnames(data.filter)
group_info2$samples <- rownames(group_info2)

pca_sample <- cbind(pca_sample, group_info2)
pca_sample  #作图数据中包含了样本坐标和分组信息

#ggplot2 绘制二维散点图
library(ggplot2)
col5 <- pal_npg(palette = c("nrc"))(7)

p <- ggplot(data = pca_sample, aes(x = Dim.1, y = Dim.2)) +
					geom_point(aes(color = group_info), size = 5) +  #根据样本坐标绘制二维散点图
					scale_color_manual(values = col5) +  #自定义颜色
					theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent',linewidth = 2), legend.key = element_rect(fill = 'transparent'), #去除背景和网格线、设置边框线粗细
							text = element_text (size = 18),axis.text = element_text (size = 18)) +  #设置字体大小
					labs(x =  paste('PCA1:', pca_eig1, '%'), y = paste('PCA2:', pca_eig2, '%'), color = '')  #将 PCA 轴贡献度添加到坐标轴标题中

p

library(plyr)
cluster_border <- ddply(pca_sample, 'group_info', function(df) df[chull(df[[1]], df[[2]]), ])
p2 <- p + geom_polygon(data = cluster_border, aes(fill = group_info), alpha = 0.2, show.legend = FALSE) +
scale_fill_manual(values = col5)
p2
ggsave("00_BatchEffect_PCA_orignial.pdf", plot = p2, width = 8, height = 5)


#3. sva包抹除批次效应
# https://cloud.tencent.com/developer/article/2319682  有ComBat和ComBat_seq可选，前者更适合芯片数据，后者针对测序数据

#sva包ComBat for MicroArray data
#library(sva)
#Expr <- data.filter
#coldata = group_info2
#design <- model.matrix(~coldata$group_info)
#combat_Expr <-ComBat(dat = Expr,batch = coldata$batch)

#sva包 ComBat_seq for NGS data
#计算每个样品的counts总数
library(sva)
a <- group_info2
a$sumcounts = 0
for (i in 1:length(colnames(data.filter)))
{
  a[i,"sumcounts"] = sum(data.filter[,i])
}
write.csv(a, file="00_data_for_batch_identification.csv")
##导出文件到Excel，根据counts数据大小设置batch，然后保存

b <- read.csv("00_data_for_batch_identification.csv", header = T,row.names = 1)
coldata = b
coldata
combat_Expr <- ComBat_seq(counts = as.matrix(data.filter), 
                                batch = coldata$batch,
								group = coldata$group_info
)

write.csv(combat_Expr, file="newdata_filter_remove_pre_batheffect_removed.csv")

#相关性分析
#计算距离
sample_cor <- cor(combat_Expr)
sample_cor1 <- round(sample_cor, digits = 2)
#画图
library(pheatmap)
cor_heatmap <- pheatmap(sample_cor1, display_numbers = T,fontsize = 10, angle_col = 45)
cor_heatmap
ggsave("00_BatchEffect_Correlation_sva_heatmap.pdf", plot = cor_heatmap, width = 11, height = 8) 


#PCA 分析
gene <- t(combat_Expr)
library(FactoMineR)
gene.pca2 <- PCA(gene, ncp = 2, scale.unit = TRUE, graph = FALSE)
plot(gene.pca)  #PCA 简图
#提取样本在 PCA 前两轴中的坐标
pca_sample2 <- data.frame(gene.pca2$ind$coord[ ,1:2])
head(pca_sample2)
pca_sample2
#提取 PCA 前两轴的贡献度
pca_eig1 <- round(gene.pca2$eig[1,2], 2)
pca_eig2 <- round(gene.pca2$eig[2,2],2 )
pca_eig1
pca_eig2

#增添分组信息
colnames(combat_Expr)
group_info = c(rep("Ras_GFP_neg",3), rep("Ras_GFP_pos",3), rep("RasPPP1CA_GFP_neg",3), 
rep("RasPPP1CA_GFP_pos",4), rep("WT_GFP_neg",3), rep("WT_GFP_pos",3))
group_info2 <- as.data.frame(group_info)
rownames(group_info2) <- colnames(combat_Expr)
group_info2$samples <- rownames(group_info2)
#合并pca_sample2和group_info2
pca_sample2 <- cbind(pca_sample2, group_info2)
pca_sample2  #作图数据中包含了样本坐标和分组信息
#ggplot2 绘制二维散点图
library(ggplot2)
col5 <- pal_npg(palette = c("nrc"))(7)

p <- ggplot(data = pca_sample2, aes(x = Dim.1, y = Dim.2)) +
					geom_point(aes(color = group_info), size = 5) +  #根据样本坐标绘制二维散点图
					scale_color_manual(values = col5) +  #自定义颜色
					theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent',linewidth = 2), legend.key = element_rect(fill = 'transparent'), #去除背景和网格线、设置边框线粗细
							text = element_text (size = 18),axis.text = element_text (size = 18)) +  #设置字体大小
					labs(x =  paste('PCA1:', pca_eig1, '%'), y = paste('PCA2:', pca_eig2, '%'), color = '')  #将 PCA 轴贡献度添加到坐标轴标题中

p
#多边形连接同类别对象边界的样式，适用于各组样本数大于 3 个的情况
library(plyr)
cluster_border <- ddply(pca_sample2, 'group_info', function(df) df[chull(df[[1]], df[[2]]), ])
p2 <- p + geom_polygon(data = cluster_border, aes(fill = group_info), alpha = 0.2, show.legend = FALSE) +
scale_fill_manual(values = col5)

p2
ggsave("00_BatchEffect_PCA_sva.pdf", plot = p2, width = 8, height = 5)


#读取去批次效应之后的Count矩阵
combat_Expr <- read.csv("newdata_filter_remove_pre_batheffect_removed.csv", header = T,row.names = 1)

# condition table
colnames(combat_Expr)
group_info = c(rep("Ras_GFP_neg",3), rep("Ras_GFP_pos",3), rep("RasPPP1CA_GFP_neg",3), 
rep("RasPPP1CA_GFP_pos",4), rep("WT_GFP_neg",3), rep("WT_GFP_pos",3))

#DGEList constructor
dge.list.obj <- DGEList(counts = combat_Expr, group = group_info)
dge.list.obj
# Normalization method: "TMM","TMMwsp","RLE","upperquartile","none"
dge.list.obj <- calcNormFactors(dge.list.obj,method = "RLE") # DESeq2, cuffdiff
dge.list.obj$samples
#plotMDS(dge.list.obj)
# make design matrix
design.mat <- model.matrix(~group_info)
# estimate dispersion
dge.list.obj <- estimateDisp(dge.list.obj,design.mat)
dge.list.obj$common.dispersion
dge.list.obj$tagwise.dispersion
	### Equal to the following steps
	#### 1st common dispersion
	#dge.list.obj <- estimateCommonDisp(dge.list.obj)
	#### 2nd tagwise dispersion
	#dge.list.obj <- estimateTagwiseDisp(dge.list.obj)
	#### plot dispersion
	#plotBCV(dge.list.obj, cex = 0.8)
	#### plot var and mean
	#plotMeanVar(dge.list.obj, show.raw=TRUE, show.tagwise=TRUE, show.binned=TRUE)

# test with likelihood ratio test
#1.  GLM general linear model
	#fit <- glmFit(dge.list.obj, design.mat)
	#lrt <- glmLRT(fit, coef=2)
	#DEGs.res.lrt <- as.data.frame(topTags(lrt,n=nrow(count_df.filter),sort.by = "logFC"))
#2.  exactTest
dge.list.res <- exactTest(dge.list.obj, pair = c("WT_GFP_neg","Ras_GFP_neg")) #control在前
#Down None   Up 
#1361 8205 1481

dge.list.res <- exactTest(dge.list.obj, pair = c("Ras_GFP_neg","RasPPP1CA_GFP_neg")) #control在前
#Down None   Up 
#1324 8779	944

dge.list.res <- exactTest(dge.list.obj, pair = c("WT_GFP_neg","RasPPP1CA_GFP_neg")) #control在前
#Down None   Up 
# 907 9582  558

dge.list.res <- exactTest(dge.list.obj, pair = c("WT_GFP_pos","Ras_GFP_pos")) #control在前
#Down None   Up 
#2681 4791 3575

dge.list.res <- exactTest(dge.list.obj, pair = c("Ras_GFP_pos","RasPPP1CA_GFP_pos")) #control在前
#Down None   Up 
#3090 5267 2690

dge.list.res <- exactTest(dge.list.obj, pair = c("WT_GFP_pos","RasPPP1CA_GFP_pos")) #control在前
#Down None   Up 
#1080 7813 2154

#topTags 按不同的标准来排列差异数据
topTags(dge.list.res)
ordered_tags <- topTags(dge.list.res, n=1000000, sort.by = "logFC")
head(ordered_tags)

DEGs.res <- as.data.frame(topTags(dge.list.res,n=nrow(combat_Expr),sort.by = "logFC"))

# 输出相关的文件
allDiff=DEGs.res[is.na(DEGs.res$FDR)==FALSE,]
diff=allDiff
###对基因进行注释-获取gene_symbol
library("biomaRt")
library("curl")
mart <- useDataset("dmelanogaster_gene_ensembl", useMart("ensembl"))
my_ensembl_gene_id<-row.names(diff)
study_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"), filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
head(study_symbols)
ensembl_gene_id<-rownames(diff)
diff <- cbind(ensembl_gene_id,diff)
colnames(diff)[1]<-c("ensembl_gene_id")
#全部基因的edger输出结果--diff_name
diff_name <- merge(diff,study_symbols,by="ensembl_gene_id")
write.table(diff_name,file="Part1_edgerOut.xls",sep="\t",quote=F)
#所有的差异基因--diff_nameSig
diff_nameSig = diff_name[(diff_name$FDR < padj & (diff_name$logFC>foldChange | diff_name$logFC<(-foldChange))),]
write.table(diff_nameSig, file="Part1_diff_nameSig.xls",sep="\t",quote=F)
nrow(diff_nameSig)
#所有的上调基因--diff_nameUp
diff_nameUp = diff_name[(diff_name$FDR < padj & (diff_name$logFC>foldChange)),]
write.table(diff_nameUp, file="Part1_diff_nameup.xls",sep="\t",quote=F)
nrow(diff_nameUp)
#所有的下调基因--diff_nameDown
diff_nameDown = diff_name[(diff_name$FDR < padj & (diff_name$logFC<(-foldChange))),]
write.table(diff_nameDown, file="Part1_diff_namedown.xls",sep="\t",quote=F)
nrow(diff_nameDown)
#所有基因normalized counts table--normalizeExp
normalizeExp=cpm(dge.list.obj)
write.table(normalizeExp,file="Part1_normalizeExp.xls",sep="\t",quote=F,col.names=F)   #输出所有基因校正后的表达值（normalizeExp.txt）

#DEGs normalized counts table--diffExp
diffExp=normalizeExp[diff_nameSig$ensembl_gene_id,]
write.table(diffExp,file="Part1_diffmRNAExp.xls",sep="\t",quote=F,col.names=F)         #输出差异基因校正后的表达值（diffmRNAExp.txt）
nrow(diffExp)

# MA plot
select.sign.gene_id = rownames(dge.list.res)[as.logical(rownames(diffExp))]

plotSmear(dge.list.res, de.tags = rownames(diffExp), cex = 0.5,ylim=c(-5,5)) 

# Step 2 Volcano Plot
#画火山图
library("ggplot2")
library("ggpubr")
library("ggthemes")
library("pheatmap")
#在画图之前需要将padj转换成-1*log10，这样可以拉开表达基因之间的差距    -log10(FDR=0.05)=1.30103 log2(FC=1.5)=0.5849625 -log10(FDR=0.01)=2
diff_name2 <- diff_name
rownames(diff_name2)=diff_name2[,1]
diff_name2$log10diff_namepadj <- -log10(diff_name2$FDR)
diff_name2$group <- "nonsignificance"
diff_name2$group[(diff_name2$log10diff_namepadj > 1.30103) & (diff_name2$logFC > 0.5849625)]="Up"
diff_name2$group[(diff_name2$log10diff_namepadj > 1.30103) & (diff_name2$logFC < -0.5849625)]="Down"
table(diff_name2$group)
#Label-FDR最高的20个基因(ensembl_gene_id和external_gene_name可以互换)
diff_name2$label=""
diff_name2 <- diff_name2[order(diff_name2$FDR),]
diff_name2_upgenes_FDR <- head(diff_name2$external_gene_name[which(diff_name2$group=="Up")],20)
diff_name2_downgenes_FDR <- head(diff_name2$external_gene_name[which(diff_name2$group=="Down")],20)
diff_name2_top20genes_FDR <- c(as.character(diff_name2_upgenes_FDR),as.character(diff_name2_downgenes_FDR))
diff_name2$label[match(diff_name2_top20genes_FDR,diff_name2$external_gene_name)] <- diff_name2_top20genes_FDR
#Label-logFC最高的10个基因
diff_name2$logFC_abs=abs(diff_name2$logFC)
diff_name2 <- diff_name2[order(diff_name2$logFC_abs),]
diff_name2_upgenes_logFC <- tail(diff_name2$external_gene_name[which(diff_name2$group=="Up")],10)
diff_name2_downgenes_logFC <- head(diff_name2$external_gene_name[which(diff_name2$group=="Down")],10)
diff_name2_top10genes_logFC <- c(as.character(diff_name2_upgenes_logFC),as.character(diff_name2_downgenes_logFC))
diff_name2$label[match(diff_name2_top10genes_logFC, diff_name2$external_gene_name)] <- diff_name2_top10genes_logFC
write.table(diff_name2,file="Part1_diff_name2.xls",sep="\t",quote=F,col.names=T)   
#画分界线（1.30是-log10(FDR=0.05)的值） log10(FDR=0.05)=-1.30103 log2(FC=1.5)=0.5849625
volcano <- ggscatter(data = diff_name2,x = "logFC",y = "log10diff_namepadj",
	color = "group", 
	palette = c("#2f5688","#BBBBBB","#CC0000"),
	size = 1,
	font.label = c(8, "plain"),
	label = diff_name2$label,
	repel = T,
	xlab="Log2FoldChange",
	ylab="-Log10(FDR)",)+theme_base()+
	geom_hline(yintercept = 1.3,linetype="dashed")+ 
	geom_vline(xintercept = c(-0.5849625,0.5849625),linetype="dashed")
ggsave("Part2_Volcanoplot.pdf", plot = volcano, width = 11, height = 8) 

#MAplot
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(DESeq2)
library(ggrepel)

ggplot(diff_name2,aes(x=logCPM,y=logFC)) +
    geom_point(aes(color=group, alpha = group),size = 1.2,
               show.legend = T) +
    geom_hline(yintercept =  c(-0.5849625, 0.5849625),lty=2,lwd = 1) +
    theme_bw(base_size = 12) + 
    ggtitle("Neg vs Pos") + 
    scale_color_manual(values = c("cadetblue","black","grey","darkorange")) +
#	scale_color_manual(values = c("#2f5688","grey","#CC0000","darkorange")) +
    scale_size_manual(values = c(1.5,2,1.5,1.5)) +
    scale_alpha_manual(values = c(0.6,1,0.6,0.6)) +
    scale_x_log10() + 
    theme(panel.grid=element_blank(),
	panel.border = element_rect(color = "black", size = 2, fill = NA),  #设置黑框颜色和粗细
	axis.text.x = element_text(face = "plain", size = 12, colour = "black"),   #字体调整  face = "italic" ("plain", "italic", "bold", "bold.italic")
	axis.text.y =element_text(face = "plain", size = 12, colour = "black")) +   
    scale_y_continuous(limits = c(-10,10),breaks = c(-10,-5,0,5,10)) +
    xlab("Mean Normalized Counts") + 
    ylab("Log2FoldChange") +
    geom_text_repel(aes(label=NA), color="black",fontface="italic",  #若要Label gene id则为label=label
                    size=4, segment.size=0.5,hjust=0,
                    nudge_x=1, nudge_y = 1)
ggsave(filename="Part2_MAplot.pdf", plot=last_plot(), scale = 1, width = 13, height = 8, units = c("cm"), dpi = 600,device = "pdf") 

#cadetblue   #9EC4C5
#grey   #DDDDDD


# Step 2 Volcano Plot New
#From edgeR results to count_matrix
DEGs.res <- as.data.frame(topTags(dge.list.res,n=nrow(combat_Expr),sort.by = "logFC"))
# 输出相关的文件
allDiff=DEGs.res[is.na(DEGs.res$FDR)==FALSE,]
diff=allDiff
###对基因进行注释-获取gene_symbol
library("biomaRt")
library("curl")
mart <- useDataset("dmelanogaster_gene_ensembl", useMart("ensembl"))
my_ensembl_gene_id<-row.names(diff)
study_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"), filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
head(study_symbols)
ensembl_gene_id<-rownames(diff)
diff <- cbind(ensembl_gene_id,diff)
colnames(diff)[1]<-c("ensembl_gene_id")
#全部基因的edger输出结果--diff_name
diff_name <- merge(diff,study_symbols,by="ensembl_gene_id")

#准备需要的文件diff_name2--input
diff_name2 <- diff_name
rownames(diff_name2)=diff_name2[,1]
diff_name2$log10diff_namepadj <- -log10(diff_name2$FDR)
diff_name2$group <- "None"
diff_name2$group[(diff_name2$log10diff_namepadj > 1.30103) & (diff_name2$logFC > 0.5849625)]="Up"
diff_name2$group[(diff_name2$log10diff_namepadj > 1.30103) & (diff_name2$logFC < -0.5849625)]="Down"
table(diff_name2$group)

#需要label的gene
genes = c("Dif","dl","Myd88","NT1","pll","spz","spz5","Tehao","Tl","Toll-7","tub","wek","Def","Mtk","Dro","Drs","Bombc2")
genes = c("AttD","BomBc1","BomBc2","BomS6","BomS5","BomS2","BomS1","BomS3","Mtk","DptB","DptA","Def","Dro","Drs")
genes = c("Pvf1","Pvf2","Ets21C")
genes = c("")
#label genes of interest
diff_name2$label = ifelse(diff_name2$external_gene_name %in% genes, diff_name2$external_gene_name,"")

#diff_name2 to input
head(diff_name2)
input <- diff_name2[,c("ensembl_gene_id","external_gene_name","logFC","FDR","log10diff_namepadj","group","label")]
table(input$group)

#如何让你的火山图标签“乖乖听话”？个性化标签全拿捏！https://mp.weixin.qq.com/s/JpA7egexar9oN0QFjVt2EQ
#自定义颜色：
mycol <- c("#00468B","#d8d8d8","#EB4232")
#自定义主题：
mytheme <- theme_classic() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        plot.margin = margin(15,5.5,5.5,5.5))

#ggplot2绘制火山图：
##foldChange=0.5849625    foldChange=log2(FC=1.5)=0.5849625
foldChange=0.5849625
padj=0.05
p <- ggplot(data = input,
            aes(x = logFC,
                y = log10diff_namepadj,
                color = group)) + #建立映射
  geom_point(size = 2.2) + #绘制散点
  scale_colour_manual(name = "", values = alpha(mycol, 0.7)) + #自定义散点颜色
  scale_x_continuous(limits = c(-15, 15),
                     breaks = seq(-15, 15, by = 3)) + #x轴限制
  scale_y_continuous(expand = expansion(add = c(2, 0)),
                     limits = c(0, 50),
                     breaks = seq(0, 50, by = 5)) + #y轴限制
  xlab("Log2FoldChange") + #X轴标签
  ylab("-Log10FDR") + #Y轴标签
  geom_hline(yintercept = c(-log10(padj)),size = 0.7,color = "black",lty = "dashed") + #水平阈值线
  geom_vline(xintercept = c(-foldChange, foldChange),size = 0.7,color = "black",lty = "dashed") + #垂直阈值线
  mytheme
p

labeled_id <- input[(input$label != ""),]
labeled_id
head(labeled_id)

p2 <- p +
    geom_point(data = labeled_id,
               aes(x = logFC, y = log10diff_namepadj),
               color = '#EB4232', size = 4.5, alpha = 0.2) + #图中强调显示一下需要标注的散点
    geom_text_repel(data = labeled_id,
                    aes(x = logFC, y = log10diff_namepadj, label = label),
                    seed = 233, #可设置随机数种子
                    size = 3, #字体大小
                    color = 'black', #字体颜色
                    min.segment.length = 0, #始终为标签添加指引线段；若不想添加线段，则改为Inf
                    force = 2, #重叠标签间的排斥力
                    force_pull = 2, #标签和数据点间的吸引力
                    box.padding = 0.4, #标签周边填充量，默认单位为行
                    max.overlaps = Inf, #排斥重叠过多标签，设置为Inf则可以保持始终显示所有标签
                    segment.linetype = 3, #线段类型,1为实线,2-6为不同类型虚线
                    segment.color = 'black', #线段颜色
                    segment.alpha = 0.8, #线段不透明度
                    nudge_x = 15 - labeled_id$logFC, #标签x轴起始位置调整
                    nudge_y = 3 + labeled_id$log10diff_namepadj, #标签y轴起始位置调整
                    direction = "y", #按y轴调整标签位置方向，若想水平对齐则为x
                    hjust = 1 #对齐标签：0右对齐，1左对齐，0.5居中
    )
p2

ggsave(filename = "Part2_VolcanoPlotNew_Labeled.pdf",plot = p2, width = 13,height = 12, units = "cm", dpi = 600)

#ggsave(filename = "Part2_VolcanoPlotNew_Labeled_Pvf.pdf",plot = p2, width = 15,height = 11, units = "cm", dpi = 600)


#分别筛选上下调中显著性top20，作为本次测试需要添加的目标标签（共40个标签）：
up <- filter(input, group == 'Up') %>% top_n(20, log10diff_namepadj) %>% summarise(ensembl_gene_id,
                                                                                   external_gene_name,
                                                                                   logFC,
                                                                                   FDR,
                                                                                   log10diff_namepadj,
                                                                                   group,
                                                                                   label=external_gene_name)
down <- filter(input, group == 'Down') %>% top_n(20, log10diff_namepadj) %>% summarise(ensembl_gene_id,
                                                                                       external_gene_name,
                                                                                       logFC,
                                                                                       FDR,
                                                                                       log10diff_namepadj,
                                                                                       group,
                                                                                       label=external_gene_name)
head(up);head(down)

#两个函数都可用于添加标签：
##geom_text_repel():常规样式，默认参数配置
#p1 <- p +
#  geom_text_repel(data = up,
#                  aes(x = log2FoldChange, y = -log10(pvalue), label = Symbol)) +
#  geom_text_repel(data = down,
#                  aes(x = log2FoldChange, y = -log10(pvalue), label = Symbol))
#p1
##geom_label_repel():带边框样式，默认参数配置
#p2 <- p +
#  geom_label_repel(data = up,
#                   aes(x = log2FoldChange, y = -log10(pvalue), label = Symbol)) +
#  geom_label_repel(data = down,
#                   aes(x = log2FoldChange, y = -log10(pvalue), label = Symbol))
#p2

p5 <- p +
  geom_point(data = up,
             aes(x = logFC, y = log10diff_namepadj),
             color = '#EB4232', size = 4.5, alpha = 0.2) + #图中强调显示一下需要标注的散点
  geom_text_repel(data = up,
                  aes(x = logFC, y = log10diff_namepadj, label = label),
                  seed = 233, #可设置随机数种子
                  size = 3.5, #字体大小
                  color = 'black', #字体颜色
                  min.segment.length = 0, #始终为标签添加指引线段；若不想添加线段，则改为Inf
                  force = 2, #重叠标签间的排斥力
                  force_pull = 2, #标签和数据点间的吸引力
                  box.padding = 0.4, #标签周边填充量，默认单位为行
                  max.overlaps = Inf, #排斥重叠过多标签，设置为Inf则可以保持始终显示所有标签
                  segment.linetype = 3, #线段类型,1为实线,2-6为不同类型虚线
                  segment.color = 'black', #线段颜色
                  segment.alpha = 0.8, #线段不透明度
                  nudge_x = 9 - up$logFC, #标签x轴起始位置调整
                  direction = "y", #按y轴调整标签位置方向，若想水平对齐则为x
                  hjust = 1 #对齐标签：0右对齐，1左对齐，0.5居中
  )
p5

#左侧添加方法相同：
p6 <- p5 +
  geom_point(data = down,
             aes(x = logFC, y = log10diff_namepadj),color = '#2DB2EB',size = 4.5, alpha = 0.2) +
  geom_text_repel(data = down,
                  aes(x = logFC, y = log10diff_namepadj, label = label),
                  seed = 233,
                  size = 3.5,
                  color = 'black',
                  min.segment.length = 0,
                  force = 2,
                  force_pull = 2,
                  box.padding = 0.4,
                  max.overlaps = Inf,
                  segment.linetype = 3,
                  segment.color = 'black',
                  segment.alpha = 0.8,
                  nudge_x = -8 - down$logFC,
                  direction = "y",
                  hjust = 1 #改为左对齐即可
  )
p6
ggsave(filename = "Part2_VolcanoPlotNew.pdf",plot = p6, width = 14,height = 12, units = "cm", dpi = 600)



##Step 3 DEGs Heatmap
library(devtools)
library(biomaRt)
library(curl)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(stringr)
library(edgeR)
library(ggsci)
library(cowplot)
library(tidyverse)
library(ggunchull)
library(SCENIC)
library(scales)
setwd("H:\\PPP1CA&PPP1CC\\01BulkRNAseq\\02Analysis")                    #设置工作目录

#调取标准化之后的基因表达矩阵
normalizeExp=cpm(dge.list.obj)
colnames(normalizeExp)

#Heatmap

#genes of interest
#genes = c("DptA","DptB","AttA","AttB","AttC","AttD","Dro","CecA1","CecA2","CecB","CecC","Anp","Def","Drs","Drsl1","Drsl2","Drsl3","Drsl4","Drsl5","Drsl6","Mtk","BomBc1","BomBc2","BomBc3","BomT1","BomT2","BomT3","BomS1","BomS2","BomS3","BomS4","BomS5","BomS6")
genes <- read.table(file="DEGs_GFP_neg_7+33candidates.txt",header = F)
normalizeExp_neg = normalizeExp[,c("WT_GFP_neg13","WT_GFP_neg14","WT_GFP_neg15",
															"Ras_GFP_neg1","Ras_GFP_neg2","Ras_GFP_neg3",
															"RasPPP1CA_GFP_neg1","RasPPP1CA_GFP_neg2","RasPPP1CA_GFP_neg3")]

head(normalizeExp_neg)
Log_cpm_1_normalizeExp <- as.matrix(log2(normalizeExp_neg+1))

genes <- read.table(file="DEGs_GFP_pos_112+12candidates.txt",header = F)
normalizeExp_pos = normalizeExp[,c("WT_GFP_pos16","WT_GFP_pos17","WT_GFP_pos18",
															"Ras_GFP_pos4","Ras_GFP_pos5","Ras_GFP_pos6",
															"RasPPP1CA_GFP_pos4","RasPPP1CA_GFP_pos5","RasPPP1CA_GFP_pos6","RasPPP1CA_GFP_pos7")]
head(normalizeExp_pos)
Log_cpm_1_normalizeExp <- as.matrix(log2(normalizeExp_pos+1))
###对基因进行注释-获取gene_symbol
library("biomaRt")
library("curl")
mart <- useDataset("dmelanogaster_gene_ensembl", useMart("ensembl"))
study_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"), filters = 'external_gene_name', values = genes, mart = mart)
head(study_symbols)

#Heatmap 画热图
library(pheatmap)
library(RColorBrewer)
#可以用的多种颜色
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
colors <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
colors <- colorRampPalette(c("cadetblue","black","grey"))(50)
#特定list里全部基因的热图
#调取所有基因的cpm文件

#Log_cpm_1_normalizeExp <- as.matrix(normalizeExp)
Log_cpm_1_normalizeExp <- as.data.frame(Log_cpm_1_normalizeExp)
head(Log_cpm_1_normalizeExp)

input_heatmap <- Log_cpm_1_normalizeExp[study_symbols$ensembl_gene_id,]
input_heatmap2 <- na.omit(input_heatmap)

rownames(study_symbols) <- study_symbols[,1]
rownames(input_heatmap2) <- study_symbols[rownames(input_heatmap2),"external_gene_name"] 

input_heatmap2
head(input_heatmap2)
#input_heatmap_regroup = input_heatmap2[,c("Ras1A", "Ras2A", "Ras3A", "M6-1B", "M6-2B", "M6-3B")]
input_heatmap_regroup = input_heatmap2  #这里不用调整
write.table(input_heatmap_regroup,file="Part3_input_heatmap_regroup.xls",sep="\t",quote=F,col.names=T)         #输出差异基因校正后的表达值
nrow(input_heatmap_regroup)

#竖版
colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
colors <- colorRampPalette(c("cadetblue","black","lightgrey"))(50)
heatmap=pheatmap(input_heatmap_regroup,color = colors,
                 main="",
                 fontsize = 12,
                 scale="row",
                 border_color = "black",
                 na_col = "grey",
                 cluster_rows = T,cluster_cols = F,
                 show_rownames = T,show_colnames = T,
                 treeheight_row = 20,treeheight_col = 10,
                 cellheight = 10,cellwidth = 15,
                 cutree_row=2,cutree_col=2,
                 display_numbers = F,legend = T,
)
heatmap
ggsave("Part3_Heatmap_genes_Vertical.pdf", plot = heatmap, width = 8, height = 22) 

#横版
colors <- colorRampPalette(c("cadetblue","black","lightgrey"))(50)
heatmap=pheatmap(t(input_heatmap_regroup),color = colors,
                 main="",
                 fontsize = 12,
                 scale="column",
                 border_color = "black",
                 na_col = "grey",
                 cluster_rows = F,cluster_cols = T,
                 show_rownames = T,show_colnames = T,
                 treeheight_row = 20,treeheight_col = 10,
                 cellheight = 10,cellwidth = 10,
                 cutree_row=2,cutree_col=2,
                 display_numbers = F,legend = T,
)
heatmap
ggsave("Part3_Heatmap_genes_Horizontal.pdf", plot = heatmap, width = 22, height = 8) 















# Step 4 差异表达基因功能富集分析
library(clusterProfiler)
library(DOSE)
#library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(ggplot2)
library(stringr)
library(AnnotationDbi)
library(org.Dm.eg.db)
library(Cairo)
library(enrichplot)
#GO Enrichment 
head(diff_nameSig)
#keyType = "ENSEMBL", 这里只支持ENSEMBL
ALL <- enrichGO(gene = diff_nameSig$ensembl_gene_id, 
                OrgDb = org.Dm.eg.db, 
                keyType = "ENSEMBL",
                ont = 'ALL',
                pvalueCutoff  = 0.05,
                pAdjustMethod = "BH",  
                qvalueCutoff  = 0.1, readable=T)  #一步到位
				
write.csv(as.data.frame(ALL@result), file="Part2_GOALL.csv",sep="\t", quote=FALSE)
GOFile = read.csv("Part2_GOALL.csv", header = T, sep="\t")

Goplot <- dotplot(ALL, split = "ONTOLOGY", font.size = 8, showCategory = 10) + 
				facet_grid(ONTOLOGY ~ ., scale = "free") + 
				scale_y_discrete(labels = function(x) str_wrap(x, width =50)) + 
                scale_size(range=c(2, 6))   #设置点的大小

ggsave("Part2_Goplot.pdf", plot = Goplot, width = 8, height = 8) 
head(ALL,1);dim(ALL)				
sum(ALL$ONTOLOGY=="BP") #Biological process基因产物参与的生物路径或机制
sum(ALL$ONTOLOGY=="CC") #Cellular component基因产物在细胞内外的位置
sum(ALL$ONTOLOGY=="MF") #Molecular function基因产物分子层次的功能

################################################################################error in kegg
library(stringr)
library(DOSE)
data(geneList, package="DOSE")
columns(org.Dm.eg.db)
gene.df <- bitr(diff_nameSig$ensembl_gene_id, fromType = "ENSEMBL", 
              toType = c("SYMBOL","ENTREZID"),
              OrgDb = org.Dm.eg.db) 
			  
gene.kegg <- bitr_kegg(gene.df$ENTREZID,fromType="ncbi-geneid",
                        toType="kegg",organism='dme')
head(gene.kegg)

#options(clusterProfiler.download.method = "auto")

ekegg <- enrichKEGG(gene = gene.kegg$kegg, 
                    organism = "dme", 
                    keyType = "kegg",
                    pvalueCutoff  = 0.05,
                    pAdjustMethod = "BH",  
                    qvalueCutoff  = 0.1)  #一步到位

head(ekegg@result)
write.csv(as.data.frame(ekegg@result), file="Part2_KEGGALL.csv",sep="\t", quote=FALSE)


erich.go.BP = enrichGO(gene = DEG.gene_symbol,
                       OrgDb = org.Dm.eg.db,
                       keyType = "SYMBOL",
                       ont = "BP",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)

erich.go.CC = enrichGO(gene = DEG.gene_symbol,
                       OrgDb = org.Dm.eg.db,
                       keyType = "SYMBOL",
                       ont = "CC",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)

erich.go.MF = enrichGO(gene = DEG.gene_symbol,
                       OrgDb = org.Dm.eg.db,
                       keyType = "SYMBOL",
                       ont = "MF",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)

dotplot(erich.go.CC)
dotplot(erich.go.BP)
dotplot(erich.go.MF)

# save image to file
pdf(file="./20200926-enrich.go.BP.Dotplot.pdf",width = 10,height = 6	)
dotplot(erich.go.BP)
dev.off()

# KEGG analysis
#查看数据库中的所有名字
columns(org.Dm.eg.db)

# convert id
DEG.gene_symbol$DEG.entrez_id = mapIds(x = org.Dm.eg.db,
                       keys = DEG.gene_symbol$DEG.gene_symbol,
                       keytype = "SYMBOL",
                       column = "ENTREZID")

DEG.gene_symbol$DEG.uniprot_id = mapIds(x = org.Dm.eg.db,
                       keys = DEG.gene_symbol$DEG.gene_symbol,
                       keytype = "SYMBOL",
                       column = "UNIPROT")

erich.kegg.res <- enrichKEGG(gene = DEG.gene_symbol$DEG.entrez_id,
                             organism = "dme",
                             keyType = "kegg",
                             pAdjustMethod = "BH")

barplot(erich.kegg.res)








