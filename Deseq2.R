##Installation of Deseq2
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2")
BiocManager::install("limma")
BiocManager::install("apeglm")

lilbrary(readxl)
library (DeSeq2)
library(limma)
library(apeglm)  
library(ggplot2)

setwd('/Users/weiwang/Desktop/DEGlist/Output')
readscount <- read_excel('/Users/weiwang/Desktop/DEGlist/Input/raw data_readcount_no rep2.xlsx', sheet = "gene.description")
colData <- read_excel('/Users/weiwang/Desktop/DEGlist/Input/colData_no rep2.xlsx', sheet = "Sheet1")
condition <-factor(c(rep(c("Control", "Flag22", "Pnic", "Flag22+Pnic"), 17), "Control", "Flag22", "Pnic", rep(c("Control","Flag22", "Pnic", "Flag22+Pnic"), 10)))
timepoint <- factor(c(rep(c(rep("T1", 4), rep("T2", 4), rep("T3", 4), rep("T4", 4), rep("T5", 4), rep("T6", 4), rep("T7", 4)), 2), 
                        rep("T1", 4), rep("T2", 4), rep("T3", 4), rep("T4", 3), rep("T5", 4), rep("T6", 4), rep("T7", 4),
                          rep("T1", 4), rep("T2", 4), rep("T3", 4), rep("T4", 4), rep("T5", 4), rep("T6", 4), rep("T7", 4)))

replicate <- factor(c(rep("Rep1",28 ), rep("Rep3",27), rep("Rep0",28)))
colData
head(readscount)
condition
timepoint
replicate
dds <- DESeqDataSetFromMatrix (readscount, colData, design = ~timepoint + replicate + condition)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

#PCA analysis###
#vsdata <- vst(dds, blind=FALSE) #variance stabilizing transformation
#assay(vsdata) <- limma::removeBatchEffect(assay(vsdata), vsdata$timepoint, vsdata$condition) 
#plotPCA(vsdata, intgroup = "replicate")

vsdata <- vst(dds, blind=FALSE)
assay(vsdata) <- limma::removeBatchEffect(assay(vsdata), vsdata$timepoint, vsdata$replicate) 
plotPCA(vsdata, intgroup = "condition")

#Plot of expression values
readscount_new = assay(vsdata)
head(readscount_new)
par(cex = 0.7)
n.sample=ncol(readscount)
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
par(mfrow=c(2,2))
boxplot(readscount, col = cols, main="expression value",las=2, labels = FALSE)
boxplot(readscount_new, col = cols,main="expression value",las=2, labels = FALSE)
hist(readscount_new, main = "Readscount Histgram")

# Differential expression analysis between each treatment and control
dds_norm <- DESeq(dds, minReplicatesForReplace = Inf) #Normalization; outliers is not filtered out; same result with cookscutoff
dds_norm$condition   #make sure levels are treat/untreat，otherwise needs to explicitely say it when define results
#res_Pnic <- results(dds_norm, contrast = c("condition","Pnic","Control"), cooksCutoff = FALSE) #alpha=0.05 for padj; cookCutoff (for filter out outlier)= false cuz there are lots of outlier
#res_Flag22 <- results(dds_norm, contrast = c("condition","Flag22","Control"), cooksCutoff = FALSE) #alpha=0.05 for padj; cookCutoff (for filter out outlier)= false cuz there are lots of outlier
res_Flag22_Pnic <- results(dds_norm, contrast = c("condition","Flag22+Pnic","Control"), cooksCutoff = FALSE) #alpha=0.05 for padj; cookCutoff (for filter out outlier)= false cuz there are lots of outlier
summary(res_Pnic) 
summary(res_Flag22)
summary(res_Flag22_Pnic)

#Filter out the ones with significance with padj < 0.05
sum(res_Pnic$padj<0.05, na.rm = TRUE)
sum(res_Flag22$padj<0.05, na.rm = TRUE)
sum(res_Flag22_Pnic$padj<0.05, na.rm = TRUE)

res_Pnic_data <- merge(as.data.frame(res_Pnic),
                  as.data.frame(counts(dds_norm,normalize=TRUE)),
                  by="row.names",sort=FALSE)
res_Flag22_data <- merge(as.data.frame(res_Flag22),
                       as.data.frame(counts(dds_norm,normalize=TRUE)),
                       by="row.names",sort=FALSE)
res_Flag22_Pnic_data <- merge(as.data.frame(res_Flag22_Pnic),
                         as.data.frame(counts(dds_norm,normalize=TRUE)),
                         by="row.names",sort=FALSE)

up_PnicDEG <- subset(res_Pnic_data, padj < 0.05 & log2FoldChange > 1)
down_PnicDEG <- subset(res_Pnic_data, padj < 0.05 & log2FoldChange < -1)
up_Flag22DEG <- subset(res_Flag22_data, padj < 0.05 & log2FoldChange > 1)
down_Flag22DEG <- subset(res_Flag22_data, padj < 0.05 & log2FoldChange < -1)
up_Flag22_PnicDEG <- subset(res_Flag22_Pnic_data, padj < 0.05 & log2FoldChange > 1)
down_Flag22_PnicDEG <- subset(res_Flag22_Pnic_data, padj < 0.05 & log2FoldChange < -1)


write.csv(res_Pnic_data, "all_Pnic_no rep2.csv") #Not filtered data, for volcano
write.csv(up_PnicDEG, "up_Pnic_no rep2.csv")
write.csv(down_PnicDEG, "down_Pnic_no rep2.csv")

write.csv(res_Flag22_data, "all_Flag22_no rep2.csv") #Not filtered data, for volcano
write.csv(up_Flag22DEG, "up_Flag22_no rep2.csv")
write.csv(down_Flag22DEG, "down_Flag22_no rep2.csv")

write.csv(res_Flag22_Pnic_data, "all_Flag22_Pnic_no rep2.csv") #Not filtered data, for volcano
write.csv(up_Flag22_PnicDEG, "up_Flag22_Pnic_no rep2.csv")
write.csv(down_Flag22_PnicDEG, "down_Flag22_Pnic_no rep2.csv")

par(mar=c(8,5,2,2))
boxplot(log10(assays(dds_norm)[["cooks"]]), range=0, las=2, outbg = "green",outpch = 25)

# MA plot
dds_norm <- DESeq(dds, minReplicatesForReplace = Inf) #Normalization; outliers is not filtered out; same result with cookscutoff
resultsNames(dds_norm)  #shrink makes data more compact,don't change padj，change foldchange
res_Pnic_shrink <- lfcShrink(dds_norm, coef="condition_Pnic_vs_Control", type="apeglm") 
pdf("MAplot_res_Pnic_data_no rep2.pdf", width = 6, height = 6) 
plotMA(res_Pnic_shrink, ylim=c(-10,10), alpha=0.1, main="res_Pnic_Data_MA plot: ")
dev.off()

dds_norm <- DESeq(dds, minReplicatesForReplace = Inf) #Normalization; outliers is not filtered out; same result with cookscutoff
dds_norm$condition
res_flag22_Pnic_shrink <- lfcShrink(dds_norm, coef="condition_Flag22.Pnic_vs_Control", type="apeglm") 
pdf("MAplot_res_flag22_Pnic_data_no rep2.pdf", width = 6, height = 6) 
plotMA(res_flag22_Pnic_shrink, ylim=c(-10,10), alpha=0.1, main="res_flag22_Pnic_Data_MA plot: ")
dev.off()

dds_norm <- DESeq(dds, minReplicatesForReplace = Inf) #Normalization; outliers is not filtered out; same result with cookscutoff
dds_norm$condition
res_flag22_shrink <- lfcShrink(dds_norm, coef="condition_Flag22_vs_Control", type="apeglm") 
pdf("MAplot_res_flag22_data_no rep2.pdf", width = 6, height = 6) 
plotMA(res_flag22_shrink, ylim=c(-10,10), alpha=0.1, main="res_flag22_Data_MA plot: ")
dev.off()

voldata_Flag22 <-read.csv(file = "all_Flag22_no rep2.csv",header = TRUE, row.names =1)
voldata_Pnic <-read.csv(file = "all_Pnic_no rep2.csv",header = TRUE, row.names =1)
voldata_Flag22_Pnic <-read.csv(file = "all_Flag22_Pnic_no rep2.csv",header = TRUE, row.names =1)


pdf("volcano.pdf", width = 6.13, height = 5.18)
ggplot(data=voldata_Flag22, aes(x=log2FoldChange,y= -1*log10(padj))) +
  geom_point(aes(color='significant')) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757")) + 
  labs(title="Volcano Plot_Flag22: ", x=expression(log[2](FC), y=expression(-log[10](padj)))) +
  geom_hline(yintercept=1.3,linetype=4) +  #反对数,代表0.05的线
  geom_vline(xintercept=c(-1,1),linetype=4) +
  theme_bw() + theme(panel.grid = element_blank())  #主次网格线均为空白
ggplot(data=voldata_Pnic, aes(x=log2FoldChange,y= -1*log10(padj))) +
  geom_point(aes(color='significant')) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757")) + 
  labs(title="Volcano Plot_Pnic: ", x=expression(log[2](FC), y=expression(-log[10](padj)))) +
  geom_hline(yintercept=1.3,linetype=4) +  #反对数,代表0.05的线
  geom_vline(xintercept=c(-1,1),linetype=4) +
  theme_bw() + theme(panel.grid = element_blank())  #主次网格线均为空白
ggplot(data=voldata_Flag22_Pnic, aes(x=log2FoldChange,y= -1*log10(padj))) +
  geom_point(aes(color='significant')) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757")) + 
  labs(title="Volcano Plot_Flag22+Pnic: ", x=expression(log[2](FC), y=expression(-log[10](padj)))) +
  geom_hline(yintercept=1.3,linetype=4) +  #反对数,代表0.05的线
  geom_vline(xintercept=c(-1,1),linetype=4) +
  theme_bw() + theme(panel.grid = element_blank())  #主次网格线均为空白

dev.off()


### T4
readscount <- read_excel('/Users/weiwang/Desktop/DEGlist/Input/raw data_readcount.xlsx', sheet = "T4")
colData <- read_excel('/Users/weiwang/Desktop/DEGlist/Input/colData.xlsx', sheet = "T4")
condition <-factor(c(rep(c("Control", "Flag22", "Pnic", "Flag22+Pnic"), 2), c("Control", "Flag22", "Pnic"), c("Control", "Flag22", "Pnic", "Flag22+Pnic")))
replicate <- factor(c(rep("Rep1", 4), rep("Rep3", 3), rep("Rep0", 4)))

colData
head(readscount)
condition
replicate
dds <- DESeqDataSetFromMatrix (readscount, colData, design = ~replicate+condition)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

#PCA analysis
vsdata <- vst(dds, blind=FALSE)
assay(vsdata) <- limma::removeBatchEffect(assay(vsdata), vsdata$condition) 
plotPCA(vsdata, intgroup = "replicate")
vsdata <- vst(dds, blind=FALSE)
assay(vsdata) <- limma::removeBatchEffect(assay(vsdata), vsdata$replicate) 
plotPCA(vsdata, intgroup = "condition")

#Plot of expression values
readscount_new = assay(vsdata)
par(cex = 0.7)
n.sample=ncol(readscount)
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
par(mfrow=c(2,2))
boxplot(readscount, col = cols, main="expression value",las=2)
boxplot(readscount_new, col = cols,main="expression value",las=2)
hist(readscount_new, main = "Readscount")

# Differentially expression analysis between each treatment and control
dds_norm <- DESeq(dds, minReplicatesForReplace = Inf) #Normalization; outliers are not filtered out; same result with cookscutoff
dds_norm$condition 
res_Pnic <- results(dds_norm, contrast = c("condition","Pnic","Control"), cooksCutoff = FALSE) 
res_Flag22 <- results(dds_norm, contrast = c("condition","Flag22","Control"), cooksCutoff = FALSE) 
res_Flag22_Pnic <- results(dds_norm, contrast = c("condition","Flag22+Pnic","Control"), cooksCutoff = FALSE)

summary(res_Pnic) 
summary(res_Flag22)
summary(res_Flag22_Pnic)

#Filter out the ones with significance with padj < 0.05
sum(res_Pnic$padj<0.05, na.rm = TRUE)
sum(res_Flag22$padj<0.05, na.rm = TRUE)
sum(res_Flag22_Pnic$padj<0.05, na.rm = TRUE)

res_Pnic_data <- merge(as.data.frame(res_Pnic),
                       as.data.frame(counts(dds_norm,normalize=TRUE)),
                       by="row.names",sort=FALSE)
res_Flag22_data <- merge(as.data.frame(res_Flag22),
                         as.data.frame(counts(dds_norm,normalize=TRUE)),
                         by="row.names",sort=FALSE)
res_Flag22_Pnic_data <- merge(as.data.frame(res_Flag22_Pnic),
                              as.data.frame(counts(dds_norm,normalize=TRUE)),
                              by="row.names",sort=FALSE)

up_PnicDEG <- subset(res_Pnic_data, padj < 0.05 & log2FoldChange > 1)
down_PnicDEG <- subset(res_Pnic_data, padj < 0.05 & log2FoldChange < -1)
up_Flag22DEG <- subset(res_Flag22_data, padj < 0.05 & log2FoldChange > 1)
down_Flag22DEG <- subset(res_Flag22_data, padj < 0.05 & log2FoldChange < -1)
up_Flag22_PnicDEG <- subset(res_Flag22_Pnic_data, padj < 0.05 & log2FoldChange > 1)
down_Flag22_PnicDEG <- subset(res_Flag22_Pnic_data, padj < 0.05 & log2FoldChange < -1)


write.csv(res_Pnic_data, "T4_all_Pnic.csv") 
write.csv(up_PnicDEG, "T4_up_Pnic.csv")
write.csv(down_PnicDEG, "T4_down_Pnic.csv")

write.csv(res_Flag22_data, "T4_all_Flag22.csv") 
write.csv(up_Flag22DEG, "T4_up_Flag22.csv")
write.csv(down_Flag22DEG, "T4_down_Flag22.csv")

write.csv(res_Flag22_Pnic_data, "T4_all_Flag22_Pnic.csv") 
write.csv(up_Flag22_PnicDEG, "T4_up_Flag22_Pnic.csv")
write.csv(down_Flag22_PnicDEG, "T4_down_Flag22_Pnic.csv")

par(mar=c(8,5,2,2))
boxplot(log10(assays(dds_norm)[["cooks"]]), range=0, las=2)

resultsNames(dds_norm) 
res_shrink <- lfcShrink(dds_norm, coef=5, type="apeglm") 
pdf("T4_MAplot.pdf", width = 6, height = 6) 
plotMA(res_shrink, ylim=c(-10,10), alpha=0.1, main="T4_MA plot: ")
dev.off()

voldata_Flag22 <-read.csv(file = "T4_all_Flag22.csv",header = TRUE, row.names =1)
voldata_Pnic <-read.csv(file = "T4_all_Pnic.csv",header = TRUE, row.names =1)
voldata_Flag22_Pnic <-read.csv(file = "T4_all_Flag22_Pnic.csv",header = TRUE, row.names =1)


pdf("T4_volcano.pdf", width = 6.13, height = 5.18)
ggplot(data=voldata_Flag22, aes(x=log2FoldChange,y= -1*log10(padj))) +
  geom_point(aes(color= "significant")) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757")) + 
  labs(title="T4_Volcano Plot_Flag22: ", x=expression(log[2](FC), y=expression(-log[10](padj)))) +
  geom_hline(yintercept=1.3,linetype=4) +  #反对数,代表0.05的线
  geom_vline(xintercept=c(-1,1),linetype=4) +
  theme_bw() + theme(panel.grid = element_blank())  #主次网格线均为空白
ggplot(data=voldata_Pnic, aes(x=log2FoldChange,y= -1*log10(padj))) +
  geom_point(aes(color='significant')) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757")) + 
  labs(title="T4_Volcano Plot_Pnic: ", x=expression(log[2](FC), y=expression(-log[10](padj)))) +
  geom_hline(yintercept=1.3,linetype=4) +  #反对数,代表0.05的线
  geom_vline(xintercept=c(-1,1),linetype=4) +
  theme_bw() + theme(panel.grid = element_blank())  #主次网格线均为空白
ggplot(data=voldata_Flag22_Pnic, aes(x=log2FoldChange,y= -1*log10(padj))) +
  geom_point(aes(color='significant')) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757")) + 
  labs(title="T4_Volcano Plot_Flag22+Pnic: ", x=expression(log[2](FC), y=expression(-log[10](padj)))) +
  geom_hline(yintercept=1.3,linetype=4) +  #反对数,代表0.05的线
  geom_vline(xintercept=c(-1,1),linetype=4) +
  theme_bw() + theme(panel.grid = element_blank())  #主次网格线均为空白

dev.off()

