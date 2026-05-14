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
library(xlsx)
setwd("H:\\PPP1CA&PPP1CC\\01BulkRNAseq\\01Raw_count\\GFP")                    #设置工作目录

##foldChange=0.5849625    foldChange=log2(FC=1.5)=0.5849625
foldChange=1  
padj=0.01

#Step 1 
#整理RNA-seq的基因表达矩阵
#https://mp.weixin.qq.com/s/qsBOPoJwItO5qAGstKRC_A

#创建一个文件夹，这个folder里只有sample的count文件
#显示目前工作目录有sample文件

all_files <- list.files(path = )
print(all_files)

#提取sampleID
filesID <- as.data.frame(all_files)
print(filesID)
sampleID <- separate(data = filesID,col = all_files, into = c("id","other"),sep = "_")
print(sampleID)

#读取并合并不同sample的基因表达矩阵
#显示所有被读取的文件的路径
for (i in 1:length(all_files))
{
  files <- all_files[i]
  model_path = paste('H:\\PPP1CA&PPP1CC\\01BulkRNAseq\\01Raw_count\\GFP\\',files, sep='')
  print(model_path)
}

#将所有的sample储存在一个list对象里
q <- list()
for (i in 1:length(all_files))
{
  files <- all_files[i]
  model_path = paste('H:\\PPP1CA&PPP1CC\\01BulkRNAseq\\01Raw_count\\GFP\\',files, sep='')
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


######################
setwd("H:\\PPP1CA&PPP1CC\\01BulkRNAseq\\01Raw_count\\NonGFP")             #设置工作目录

##foldChange=0.5849625    foldChange=log2(FC=1.5)=0.5849625
foldChange=1  
padj=0.01

#Step 1 
#整理RNA-seq的基因表达矩阵
#https://mp.weixin.qq.com/s/qsBOPoJwItO5qAGstKRC_A

#创建一个文件夹，这个folder里只有sample的count文件
#显示目前工作目录有sample文件

all_files <- list.files(path = )
print(all_files)

#提取sampleID
filesID <- as.data.frame(all_files)
print(filesID)
sampleID <- separate(data = filesID,col = all_files, into = c("id","other"),sep = "_")
print(sampleID)

#读取并合并不同sample的基因表达矩阵
#显示所有被读取的文件的路径
for (i in 1:length(all_files))
{
  files <- all_files[i]
  model_path = paste('H:\\PPP1CA&PPP1CC\\01BulkRNAseq\\01Raw_count\\NonGFP\\',files, sep='')
  print(model_path)
}

#将所有的sample储存在一个list对象里
q <- list()
for (i in 1:length(all_files))
{
  files <- all_files[i]
  model_path = paste('H:\\PPP1CA&PPP1CC\\01BulkRNAseq\\01Raw_count\\NonGFP\\',files, sep='')
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

