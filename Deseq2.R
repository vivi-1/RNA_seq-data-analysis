##Installation of Deseq2
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2")

library (DeSeq2)
setwd('/Volumes/WD1/Desktop/laboratory files/Results/Altria project/DEGlist')
colData<-read.table('raw data_readcount.xlsx', header=TRUE, row.names = 1)
condition<-factor(c(rep("Mock", 7)), c(rep("Flg", 7)), c(rep("Path", 7)))
individual <- factor(c(rep("T1",3), rep("T2",3),rep("T3",3), rep("T4",3),rep("T5",3),rep("T6",4), rep("T7",4)))