library("clusterProfiler")
library(enrichR)
library(ggplot2)

#import annotation
contents <- read.csv("/Users/weiwang/Desktop/KEGG enrichment/Input/Pnic_up_KEGG.csv",header=T,sep=",")
term2name <- read.csv("/Users/weiwang/Desktop/KEGG enrichment/Input/Nb101_ID_KO.csv",header=T,sep=",")
term2gene <- read.csv("/Users/weiwang/Desktop/KEGG enrichment/Input/KO_Gene_ID.csv",header=T,sep=",")
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
barplot(xPnicUpgene,showCategory=35, includeAll=TRUE, cex.names=0.5, font.size = 6)
# bubble plot
dotplot(xPnicUpgene,showCategory=35, font.size = 6)
#enrich plot
xPnicUpgene<-enrichplot::pairwise_termsim(xPnicUpgene)
emapplot(xPnicUpgene,  showCategory = 30, cex_label_category=0.4)
g<-cnetplot(xPnicUpgene, showCategory = 20,cex_label_category = 1, cex_label_gene=0.8)
g + scale_color_manual(values = c("red","blue"))
p1 <- heatplot(xPnicUpgene)
upset(xPnicUpgene)

GOI <- read.csv("/Users/weiwang/Desktop/GoEnrichment/Output/Pnic_up_Enricher.csv",header=T,sep=",")
dotplot(as.vector(t(GOI)))

browseKEGG(xPnicUpgene, "map01100") #Metabolic pathways
browseKEGG(xPnicUpgene, "M00176")
