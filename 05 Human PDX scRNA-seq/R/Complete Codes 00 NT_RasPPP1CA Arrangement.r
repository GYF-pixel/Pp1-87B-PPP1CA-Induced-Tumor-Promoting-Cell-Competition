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
setwd("I:\\PPP1CA&PPP1CC\\03scRNA-seq_Mice\\01Analysis\\4th_Try")



#RasPPP1CA
CA1 <- Read10X(data.dir = "I:\\PPP1CA&PPP1CC\\03scRNA-seq_Mice\\01Analysis\\4th_Try\\CA1_filtered_feature_bc_matrix")
pbmc_CA1 <- CreateSeuratObject(counts = CA1,
                                min.cells = 50, 
                                min.features = 1500)
CA2 <- Read10X(data.dir = "I:\\PPP1CA&PPP1CC\\03scRNA-seq_Mice\\01Analysis\\4th_Try\\CA2_filtered_feature_bc_matrix")
pbmc_CA2 <- CreateSeuratObject(counts = CA2,
                                min.cells = 50, 
                                min.features = 1500)
CA3 <- Read10X(data.dir = "I:\\PPP1CA&PPP1CC\\03scRNA-seq_Mice\\01Analysis\\4th_Try\\CA3_filtered_feature_bc_matrix")
pbmc_CA3 <- CreateSeuratObject(counts = CA3,
                                min.cells = 50, 
                                min.features = 1500)
											
#根据每个pbmc中有多少个细胞进行标记
pbmc_CA1@meta.data$group <- rep("CA1",4810)
pbmc_CA2@meta.data$group <- rep("CA2",3467)
pbmc_CA3@meta.data$group <- rep("CA3",2364)
pbmc = merge(pbmc_CA1, y=c(pbmc_CA2, pbmc_CA3), 
             add.cell.ids = c("CA1","CA2", "CA3"),
             merge.data = TRUE)
head(pbmc@meta.data)
pbmc@meta.data$sec_group = c("CA")
table(pbmc@meta.data$sec_group)
table(pbmc@meta.data$group)

##Human Mitochondrial gene 占比
scRNAsub = pbmc						
scRNAsub[['percent.mt']] <- PercentageFeatureSet(scRNAsub,pattern = "^MT-")

#Human Ribosomal Protein是RPS和RPL开头
scRNAsub[["percent.rp"]] = PercentageFeatureSet(scRNAsub, pattern = "^RP[SL]")

#scRNAsub[["percent.rps"]] <- PercentageFeatureSet(scRNAsub, pattern = c("^RPS"))		
#scRNAsub[["percent.rpl"]] <- PercentageFeatureSet(scRNAsub, pattern = c("^RPL"))
#scRNAsub[["percent.rp"]] = scRNAsub[["percent.rps"]] + scRNAsub[["percent.rpl"]]

 
#计算红细胞比例
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNAsub@assays$RNA)) 
HB.genes <- rownames(scRNAsub@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNAsub[["percent.HB"]]<-PercentageFeatureSet(scRNAsub, features=HB.genes) 
 
violin <- VlnPlot(scRNAsub,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), 
                  group.by = "group",
                  pt.size = 0, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("vlnplot_before_qc.pdf", plot = violin, width = 12, height = 6)


scRNAsub <- subset(scRNAsub, subset = nCount_RNA > 200 & nCount_RNA < 35000 & nFeature_RNA < 6000 & nFeature_RNA > 250 & percent.mt < 20 & percent.rp < 30)


violin <- VlnPlot(scRNAsub,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), 
                  group.by = "group",
                  pt.size = 0, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("vlnplot_after_qc.pdf", plot = violin, width = 12, height = 6)

##CellCycleScoring
g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(scRNAsub))
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(scRNAsub))
scRNAsub <- CellCycleScoring(object=scRNAsub,  g2m.features=g2m_genes,  s.features=s_genes)
table(scRNAsub@meta.data[["group"]],scRNAsub@meta.data[["Phase"]])

##Split for doublets analysis
### 拆分为seurat子对象(分开pca1和pca2)
table(scRNAsub@meta.data$group)
table(scRNAsub@meta.data$sec_group)
scRNAsub.list <- SplitObject(scRNAsub, split.by = "group")
scRNAsub.list

for (i in 1:length(scRNAsub.list)) {
	scRNAsub.list[[i]] <- NormalizeData(scRNAsub.list[[i]], verbose = FALSE, normalization.method = "LogNormalize", scale.factor = 1e4)
	scRNAsub.list[[i]] <- FindVariableFeatures(scRNAsub.list[[i]], selection.method = "vst", nfeatures = 2500)
	scRNAsub.list[[i]] <- ScaleData(scRNAsub.list[[i]], vars.to.regress = c("S.Score", "G2M.Score","percent.mt", "nCount_RNA"), verbose = TRUE)
	scRNAsub.list[[i]] <- RunPCA(scRNAsub.list[[i]], features = VariableFeatures(scRNAsub.list[[i]]), npcs = 40, nfeature.print = 10, ndims.print = 1:5, verbose = T)
	pc.num=1:40
	scRNAsub.list[[i]] <- RunUMAP(scRNAsub.list[[i]], dims=pc.num)
	scRNAsub.list[[i]] <- FindNeighbors(scRNAsub.list[[i]], dims = pc.num)
	scRNAsub.list[[i]] = FindClusters(scRNAsub.list[[i]],resolution = 0.3)
###### 检测doublets 
#Doublets被定义为在相同细胞barcode下测序的两个细胞（例如被捕获在同一液滴中）
### define the expected number of doublet cellscells.
nExp <- round(ncol(scRNAsub.list[[i]]) * 0.05)  ### expect 20% doublets
### remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
#这个DoubletFinder包的输入是经过预处理（包括归一化、降维，但不一定要聚类）的 Seurat 对象
scRNAsub.list[[i]] <- doubletFinder_v3(scRNAsub.list[[i]], pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:40)	
###### 找Doublets  
### DF的名字是不固定的，因此从scRNAsub@meta.data列名中提取比较保险
DF.name = colnames(scRNAsub.list[[i]]@meta.data)[grepl("DF.classification", colnames(scRNAsub.list[[i]]@meta.data))]
###### 过滤doublet
scRNAsub.list[[i]]=scRNAsub.list[[i]][, scRNAsub.list[[i]]@meta.data[, DF.name] == "Singlet"]
#过滤到此结束
}

scRNAsub.list
reference.list <- scRNAsub.list[c("CA1","CA2", "CA3")]
scRNAsub.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:60)
scRNAsub.integrated <- IntegrateData(anchorset = scRNAsub.anchors, dims = 1:60)

# switch to integrated assay. The variable features of this assay are automatically set during
# IntegrateData
DefaultAssay(scRNAsub.integrated) <- "integrated"

# 运行标准流程并进行可视化
scRNAsub.integrated <- ScaleData(scRNAsub.integrated, verbose = FALSE)
scRNAsub.integrated <- RunPCA(scRNAsub.integrated, npcs = 60, verbose = FALSE)
scRNAsub.integrated <- RunUMAP(scRNAsub.integrated, reduction = "pca", dims = 1:60)
scRNAsub.integrated <- FindNeighbors(scRNAsub.integrated, dims = 1:60)
scRNAsub.integrated = FindClusters(scRNAsub.integrated,resolution = 0.7)

DF.name = colnames(scRNAsub.integrated@meta.data)[!grepl("pANN_0.25", colnames(scRNAsub.integrated@meta.data))]
DF.name2 = colnames(scRNAsub.integrated@meta.data)[!grepl("DF.classifications", colnames(scRNAsub.integrated@meta.data))]
DF.name3 = colnames(scRNAsub.integrated@meta.data)[!grepl("RNA_snn", colnames(scRNAsub.integrated@meta.data))]
DF.name4 = intersect(intersect(DF.name,DF.name2),DF.name3)

scRNAsub.integrated@meta.data <- scRNAsub.integrated@meta.data[,DF.name4]

saveRDS(scRNAsub.integrated, file="CCA_CA.Rds") 



#NT
NT1 <- Read10X(data.dir = "I:\\PPP1CA&PPP1CC\\03scRNA-seq_Mice\\01Analysis\\4th_Try\\NT1_filtered_feature_bc_matrix")
pbmc_NT1 <- CreateSeuratObject(counts = NT1,
                            min.cells = 50, 
                            min.features = 1500)
NT2 <- Read10X(data.dir = "I:\\PPP1CA&PPP1CC\\03scRNA-seq_Mice\\01Analysis\\4th_Try\\NT2_filtered_feature_bc_matrix")
pbmc_NT2 <- CreateSeuratObject(counts = NT2,
                            min.cells = 50, 
                            min.features = 1500)
NT3 <- Read10X(data.dir = "I:\\PPP1CA&PPP1CC\\03scRNA-seq_Mice\\01Analysis\\4th_Try\\NT3_filtered_feature_bc_matrix")
pbmc_NT3 <- CreateSeuratObject(counts = NT3,
                            min.cells = 50, 
                            min.features = 1500)
							
#根据每个pbmc中有多少个细胞进行标记
pbmc_NT1@meta.data$group <- rep("NT1",3465)
pbmc_NT2@meta.data$group <- rep("NT2",3503)
pbmc_NT3@meta.data$group <- rep("NT3",4691)

pbmc = merge(pbmc_NT1, y=c(pbmc_NT2, pbmc_NT3), 
             add.cell.ids = c("NT1","NT2", "NT3"),
             merge.data = TRUE)
head(pbmc@meta.data)
pbmc@meta.data$sec_group = c("NT")
table(pbmc@meta.data$sec_group)
table(pbmc@meta.data$group)

##Human Mitochondrial gene 占比
scRNAsub = pbmc						
scRNAsub[['percent.mt']] <- PercentageFeatureSet(scRNAsub,pattern = "^MT-")

#Human Ribosomal Protein是RPS和RPL开头
scRNAsub[["percent.rp"]] = PercentageFeatureSet(scRNAsub, pattern = "^RP[SL]")

#scRNAsub[["percent.rps"]] <- PercentageFeatureSet(scRNAsub, pattern = c("^RPS"))		
#scRNAsub[["percent.rpl"]] <- PercentageFeatureSet(scRNAsub, pattern = c("^RPL"))
#scRNAsub[["percent.rp"]] = scRNAsub[["percent.rps"]] + scRNAsub[["percent.rpl"]]

 
#计算红细胞比例
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNAsub@assays$RNA)) 
HB.genes <- rownames(scRNAsub@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNAsub[["percent.HB"]]<-PercentageFeatureSet(scRNAsub, features=HB.genes) 
 
violin <- VlnPlot(scRNAsub,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), 
                  group.by = "group",
                  pt.size = 0, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("vlnplot_before_qc.pdf", plot = violin, width = 12, height = 6)


scRNAsub <- subset(scRNAsub, subset = nCount_RNA > 200 & nCount_RNA < 35000 & nFeature_RNA < 6000 & nFeature_RNA > 250 & percent.mt < 20 & percent.rp < 30)


violin <- VlnPlot(scRNAsub,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), 
                  group.by = "group",
                  pt.size = 0, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("vlnplot_after_qc.pdf", plot = violin, width = 12, height = 6)


##CellCycleScoring


##Split for doublets analysis
### 拆分为seurat子对象(分开pca1和pca2)
table(scRNAsub@meta.data$group)
table(scRNAsub@meta.data$sec_group)
scRNAsub.list <- SplitObject(scRNAsub, split.by = "group")
scRNAsub.list

for (i in 1:length(scRNAsub.list)) {
	scRNAsub.list[[i]] <- NormalizeData(scRNAsub.list[[i]], verbose = FALSE, normalization.method = "LogNormalize", scale.factor = 1e4)
	scRNAsub.list[[i]] <- FindVariableFeatures(scRNAsub.list[[i]], selection.method = "vst", nfeatures = 2500)
	scRNAsub.list[[i]] <- ScaleData(scRNAsub.list[[i]], vars.to.regress = c("CC.Difference", "percent.mt", "nCount_RNA"), verbose = TRUE)
	scRNAsub.list[[i]] <- RunPCA(scRNAsub.list[[i]], features = VariableFeatures(scRNAsub.list[[i]]), npcs = 40, nfeature.print = 10, ndims.print = 1:5, verbose = T)
	pc.num=1:40
	scRNAsub.list[[i]] <- RunUMAP(scRNAsub.list[[i]], dims=pc.num)
	scRNAsub.list[[i]] <- FindNeighbors(scRNAsub.list[[i]], dims = pc.num)
	scRNAsub.list[[i]] = FindClusters(scRNAsub.list[[i]],resolution = 0.3)
###### 检测doublets 
#Doublets被定义为在相同细胞barcode下测序的两个细胞（例如被捕获在同一液滴中）
### define the expected number of doublet cellscells.
nExp <- round(ncol(scRNAsub.list[[i]]) * 0.05)  ### expect 20% doublets
### remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
#这个DoubletFinder包的输入是经过预处理（包括归一化、降维，但不一定要聚类）的 Seurat 对象
scRNAsub.list[[i]] <- doubletFinder_v3(scRNAsub.list[[i]], pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:40)	
###### 找Doublets  
### DF的名字是不固定的，因此从scRNAsub@meta.data列名中提取比较保险
DF.name = colnames(scRNAsub.list[[i]]@meta.data)[grepl("DF.classification", colnames(scRNAsub.list[[i]]@meta.data))]
###### 过滤doublet
scRNAsub.list[[i]]=scRNAsub.list[[i]][, scRNAsub.list[[i]]@meta.data[, DF.name] == "Singlet"]
#过滤到此结束
}

scRNAsub.list
reference.list <- scRNAsub.list[c("NT1", "NT2", "NT3")]
scRNAsub.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:60)
scRNAsub.integrated <- IntegrateData(anchorset = scRNAsub.anchors, dims = 1:60)

# switch to integrated assay. The variable features of this assay are automatically set during
# IntegrateData
DefaultAssay(scRNAsub.integrated) <- "integrated"

# 运行标准流程并进行可视化
scRNAsub.integrated <- ScaleData(scRNAsub.integrated, verbose = FALSE)
scRNAsub.integrated <- RunPCA(scRNAsub.integrated, npcs = 60, verbose = FALSE)
scRNAsub.integrated <- RunUMAP(scRNAsub.integrated, reduction = "pca", dims = 1:60)
scRNAsub.integrated <- FindNeighbors(scRNAsub.integrated, dims = 1:60)
scRNAsub.integrated = FindClusters(scRNAsub.integrated,resolution = 0.7)

DF.name = colnames(scRNAsub.integrated@meta.data)[!grepl("pANN_0.25", colnames(scRNAsub.integrated@meta.data))]
DF.name2 = colnames(scRNAsub.integrated@meta.data)[!grepl("DF.classifications", colnames(scRNAsub.integrated@meta.data))]
DF.name3 = colnames(scRNAsub.integrated@meta.data)[!grepl("RNA_snn", colnames(scRNAsub.integrated@meta.data))]
DF.name4 = intersect(intersect(DF.name,DF.name2),DF.name3)

scRNAsub.integrated@meta.data <- scRNAsub.integrated@meta.data[,DF.name4]

saveRDS(scRNAsub.integrated, file="CCA_NT_over1500.Rds") 
