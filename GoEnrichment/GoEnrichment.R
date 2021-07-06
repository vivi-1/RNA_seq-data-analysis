library("clusterProfiler")
library(enrichR)

#import annotation
contents <- read.csv("/Users/weiwang/Desktop/GoEnrichment/Pnic_up_GO.csv",header=T,sep=",")
term2name <- read.csv("/Users/weiwang/Desktop/GoEnrichment/GOterms annotations.csv",header=T,sep=",")
term2gene <- read.csv("/Users/weiwang/Desktop/GoEnrichment/GOterms_GeneList1.csv",header=T,sep=",")
term2gene <- data.frame(term2gene)
term2name <- data.frame(term2name)
head(term2name)
# Enrichment
gene <- as.factor(contents$Gene_ID)
xPnicUpgene <- enricher(gene,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.05, qvalueCutoff = 0.05)

ouf <- paste('Pnic_up_Enricher.csv',sep ="\t")
# result
write.csv(xPnicUpgene,ouf)
# barplot
barplot(xPnicUpgene,showCategory=20, includeAll=TRUE)
# bubble plot
dotplot(xPnicUpgene,showCategory=20)
