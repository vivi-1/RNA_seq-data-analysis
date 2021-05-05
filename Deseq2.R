##Installation of Deseq2
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2")
BiocManager::install("limma")
BiocManager::install("apeglm")

lilbrary("readxl")
library (DeSeq2)
library(limma)

setwd('/Volumes/WD1/Desktop/laboratory files/Results/Altria project/DEGlist')
readscount <- read_excel('/Volumes/WD1/Desktop/laboratory files/Results/Altria project/DEGlist/raw data_readcount.xlsx')
colData <- read_excel('/Volumes/WD1/Desktop/laboratory files/Results/Altria project/DEGlist/colData.xlsx')
condition <-factor(c("Mock", "Flag22","Pnic", "Flag22+Pnic"))
timepoint <- factor(c("T1","T2","T3","T4","T5","T6", "T7"))
replicate <- factor(c("One", "Two", "Three", "Four"))
colData
head(readscount)
condition
timepoint
replicate
dds <- DESeqDataSetFromMatrix (readscount, colData, design = ~timepoint + replicate + condition)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

vsdata <- vst(dds, blind=FALSE)
assay(vsdata) <- limma::removeBatchEffect(assay(vsdata), vsdata$timepoint, vsdata$replicate) 
plotPCA(vsdata, intgroup = "condition")
plotPCA(vsdata, intgroup = "replicate")

dds_norm <- DESeq(dds, minReplicatesForReplace = Inf) #标准化; 不剔除outliers; 与cookscutoff结果相同
dds_norm$condition   #保证是levels是按照后一个比前一个即trt/untrt，否则需在results时指定
res_Pnic <- results(dds_norm, contrast = c("condition","Pnic","Mock"), cooksCutoff = FALSE) #alpha=0.05可指定padj; cookCutoff是不筛选outliers因为太多了
res_Flag22 <- results(dds_norm, contrast = c("condition","Flag22","Mock"), cooksCutoff = FALSE) #alpha=0.05可指定padj; cookCutoff是不筛选outliers因为太多了
#res_Flag22_Pinic <- results(dds_norm, contrast = c("condition","Flag22+Pinic ","Mock"), cooksCutoff = FALSE) #alpha=0.05可指定padj; cookCutoff是不筛选outliers因为太多了

summary(res_Pnic) 
summary(res_Flag22)
#summary(res_Flag22_Pinic)

res_PnicOrdered <- res_Pnic[order(res_Pnic$pvalue), ] #排序
res_Flag22Ordered <- res_Flag22[order(res_Flag22$pvalue), ] #排序

sum(res_Pnic$padj<0.05, na.rm = TRUE)
sum(res_Flag22$padj<0.05, na.rm = TRUE)
#sum(res_Flag22_Pnic$padj<0.05, na.rm = TRUE)

res_Pnic_data <- merge(as.data.frame(res_Pnic),
                  as.data.frame(counts(dds_norm,normalize=TRUE)),
                  by="row.names",sort=FALSE)
res_Flag22_data <- merge(as.data.frame(res_Flag22),
                       as.data.frame(counts(dds_norm,normalize=TRUE)),
                       by="row.names",sort=FALSE)
up_PnicDEG <- subset(res_Pnic_data, padj < 0.05 & log2FoldChange > 1)
down_PnicDEG <- subset(res_Pnic_data, padj < 0.05 & log2FoldChange < -1)
up_Flag22DEG <- subset(res_Flag22_data, padj < 0.05 & log2FoldChange > 1)
down_Flag22DEG <- subset(res_Flag22_data, padj < 0.05 & log2FoldChange < -1)

write.csv(res_Pnic_data, "all_Pnic.csv") #全部基因不筛选，做火山图的背景
write.csv(up_PnicDEG, "up_Pnic.csv")
write.csv(down_PnicDEG, "down_Pnic.csv")

write.csv(res_Flag22_data, "all_Flag22.csv") #全部基因不筛选，做火山图的背景
write.csv(up_Flag22DEG, "up_Flag22.csv")
write.csv(down_Flag22DEG, "down_Flag22.csv")

par(mar=c(8,5,2,2))
boxplot(log10(assays(dds_norm)[["cooks"]]), range=0, las=2)

library(apeglm)  
resultsNames(dds_norm)  #看一下要shrink的维度;shrink数据更加紧凑,少了一项stat，但并未改变padj，但改变了foldchange
res_shrink <- lfcShrink(dds_norm, coef="condition_Pnic_vs_Mock", type="apeglm") #最推荐apeglm算法;根据resultsNames(dds)的第5个维度，coef=5，也可直接""指定;apeglm不allow contrast，所以要指定coef
pdf("MAplot.pdf", width = 6, height = 6) 
plotMA(res_shrink, ylim=c(-10,10), alpha=0.1, main="MA plot: ")
dev.off()


plotPCA(vsdata, intgroup = "condition") #自带函数