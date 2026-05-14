
#pseudobulks method
library(Seurat)
library(tidyverse)
library(dplyr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(DESeq2)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(harmony)
library(edgeR)
set.seed(1234)
getwd()
setwd("I:\\PPP1CA&PPP1CC\\03scRNA-seq_Mice\\01Analysis\\4th_Try")
setwd("H:\\PPP1CA&PPP1CC\\03scRNA-seq_Mice\\01Analysis\\4th_Try")

#数据整合
sce1 <- readRDS("CCA_CA_scRNAsub_GFP.RDS")
sce2 <- readRDS("CCA_CA_scRNAsub_mCherry.RDS")

sce1@meta.data$Group = paste(sce1@meta.data$group, "GFP", sep="_")
sce1@meta.data$sec_Group = "CA_GFP"

sce2@meta.data$Group =paste(sce2@meta.data$group, "mCherry", sep="_")
sce2@meta.data$sec_Group = "CA_mCherry"

sce3 <- readRDS("CCA_NT_scRNAsub_GFP.RDS")
sce4 <- readRDS("CCA_NT_scRNAsub_mCherry.RDS")
sce5 <- merge(x=sce3,y=sce4)
sce5@meta.data$Group = paste(sce5@meta.data$group, "GFP&mCherry", sep="_")
sce5@meta.data$sec_Group = "CA_GFP&mCherry"

sce <- merge(x=sce1,y=c(sce2,sce5))
table(sce$Group)
table(sce$sec_Group)

DefaultAssay(sce) <- "RNA"
sce <- NormalizeData(sce, verbose = FALSE, normalization.method = "LogNormalize", scale.factor = 1e4)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2500)
sce <- ScaleData(sce, vars.to.regress = c("S.Score", "G2M.Score","percent.mt", "nCount_RNA"), verbose = TRUE)
sce <- RunPCA(sce, npcs = 60, verbose = FALSE)
sce <- RunHarmony(sce, reduction.use = "pca", group.by.vars = "group", reduction.save = "harmony", verbose = FALSE)
sce <- RunUMAP(sce, reduction = "harmony", dims = 1:60)
sce <- FindNeighbors(sce, dims = 1:60)
sce = FindClusters(sce,resolution = 0.9)

saveRDS(sce, file="CCA_NT_CA_scRNAsub_GFP&mCherry.Rds") 

#样本信息, Group:  CA1_GFP CA1_mCherry     CA2_GFP CA2_mCherry     CA3_GFP CA3_mCherry
#分组信息, sec_Group: CA_GFP CA_mCherry 

bs = split(colnames(sce),sce$Group)
ct = do.call(
  cbind,lapply(names(bs), function(x){ 
    # x=names(bs)[[1]]
    kp =colnames(sce) %in% bs[[x]]
    rowSums( as.matrix(sce@assays$RNA@counts[, kp]  ))
  })
)
colnames(ct) <- names(bs)
head(ct)
ct <- as.data.frame(ct) 
ct <- select(ct,1,3,5,2,4,6:9)
head(ct)

exprSet <- ct
exprSet=exprSet[rowMeans(exprSet)>1,]
exprSet=exprSet[rowSums(exprSet)>20,]
data.filter <- exprSet
#edgeR
# condition table
colnames(data.filter)
group_info = c(rep("CA_GFP",3),rep("CA_mCherry",3),rep("NT_GFP&mCherry", 3))

#DGEList constructor
dge.list.obj <- DGEList(counts = data.filter, group = group_info)
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
fit <- glmFit(dge.list.obj, design.mat)
lrt <- glmLRT(fit, coef=2)
DEGs.res.lrt <- as.data.frame(topTags(lrt,n=nrow(data.filter),sort.by = "logFC"))
# 输出相关的文件
allDiff=DEGs.res.lrt[is.na(DEGs.res.lrt$FDR)==FALSE,]
diff=allDiff
write.csv(diff,file="scRNAseq_edgeR_LRT.csv")
#2.  exactTest
dge.list.res <- exactTest(dge.list.obj, pair = c("CA_GFP","CA_mCherry")) #第一个是Control，第二个是Case
#topTags 按不同的标准来排列差异数据
topTags(dge.list.res)
ordered_tags <- topTags(dge.list.res, n=1000000, sort.by = "logFC")
head(ordered_tags)

DEGs.res <- as.data.frame(topTags(dge.list.res,n=nrow(data.filter),sort.by = "logFC"))

# 输出相关的文件
allDiff=DEGs.res[is.na(DEGs.res$FDR)==FALSE,]
diff=allDiff
write.csv(diff,file="scRNAseq_edgeR_LRT_NT_GFP&mCherry_CA_GFP.csv")


#添加上下调信息   padj > 0.05 log2FC > log2(1.5)=0.5849625
DEG_DESeq2 <- DEG_DESeq2 %>%
  mutate(Type = if_else(padj > 0.05, 
                        "ns",
                        if_else(abs(log2FoldChange) < 0.5849625, "ns",
                                if_else(log2FoldChange >= 0.5849625, "up", "down")))) %>%
  arrange(desc(abs(log2FoldChange))) %>% rownames_to_column("Gene_Symbol")
 
 table(DEG_DESeq2$Type)
 
  # 第三步，火山图
  ggplot(DEG_DESeq2, aes(log2FoldChange,
                       -log10(padj))) +
  geom_point(size = 3.5, 
             alpha = 0.8,
             aes(color = Type),
             show.legend = T)  +
  scale_color_manual(values = c("#00468B", "gray", "#E64B35")) +
  ylim(0, 15) +
  xlim(-10, 10) +
  labs(x = "Log2(fold change)", y = "-log10(padj)") +
  geom_hline(yintercept = -log10(0.05),  
             linetype = 2,
             color = 'black',lwd=0.8) + 
  geom_vline(xintercept = c(-1, 1), 
             linetype = 2, 
             color = 'black',lwd=0.8)+theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
ggsave("volcano.pdf")










phe = unique(sce@meta.data[,c('Group','sec_Group')])#样本信息，分组信息，根据自己的metadata进行修改
phe
group_list = phe[match(names(bs),phe$Group),'sec_Group']
table(group_list)   
exprSet = ct
dim(exprSet)
exprSet=exprSet[rowMeans(exprSet)>1,]
exprSet=exprSet[rowSums(exprSet)>20,]
dim(exprSet)  
table(group_list)

#第一步，构建DEseq对象
colData <- data.frame(row.names=colnames(exprSet),group_list=group_list)
dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = colData,
                              design = ~ group_list)
							  
# 第二步，进行差异表达分析DESeq-LRT
dds2 <- DESeq(dds, 
               test = "LRT",
               reduced = ~ 1,
               fitType = "parametric",
               sfType = "poscounts",
               betaPrior = FALSE)
table(group_list)
tmp <- results(dds2,contrast=c("group_list","CA_GFP","CA_mCherry"))  #分组情况，control， case
DEG_DESeq2 <- as.data.frame(tmp[order(tmp$padj),])
head(DEG_DESeq2)
DEG_DESeq2 = na.omit(DEG_DESeq2)
head(DEG_DESeq2)
#添加上下调信息   padj > 0.05 log2FC > log2(1.5)=0.5849625
DEG_DESeq2 <- DEG_DESeq2 %>%
  mutate(Type = if_else(padj > 0.05, 
                        "ns",
                        if_else(abs(log2FoldChange) < 0.5849625, "ns",
                                if_else(log2FoldChange >= 0.5849625, "up", "down")))) %>%
  arrange(desc(abs(log2FoldChange))) %>% rownames_to_column("Gene_Symbol")
 
 table(DEG_DESeq2$Type)
 
  # 第三步，火山图
  ggplot(DEG_DESeq2, aes(log2FoldChange,
                       -log10(padj))) +
  geom_point(size = 3.5, 
             alpha = 0.8,
             aes(color = Type),
             show.legend = T)  +
  scale_color_manual(values = c("#00468B", "gray", "#E64B35")) +
  ylim(0, 15) +
  xlim(-10, 10) +
  labs(x = "Log2(fold change)", y = "-log10(padj)") +
  geom_hline(yintercept = -log10(0.05),  
             linetype = 2,
             color = 'black',lwd=0.8) + 
  geom_vline(xintercept = c(-1, 1), 
             linetype = 2, 
             color = 'black',lwd=0.8)+theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
ggsave("volcano.pdf")

