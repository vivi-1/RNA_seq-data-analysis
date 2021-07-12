#BiocManager::install("pathview")
library("clusterProfiler")
library(enrichR)
library(ggplot2)

#import annotation
contents <- read.csv("/Users/weiwang/Desktop/KEGG enrichment/Input/Pnic_up_KEGG1.csv",header=T,sep=",")
term2name <- read.csv("/Users/weiwang/Desktop/KEGG enrichment/Input/Nb101_ID_KO.csv",header=T,sep=",")
term2gene <- read.csv("/Users/weiwang/Desktop/KEGG enrichment/Input/KO_Gene_ID.csv",header=T,sep=",")
term2gene <- data.frame(term2gene)
term2name <- data.frame(term2name)
head(term2name)
# Enrichment
gene <- as.factor(contents$Gene_ID)


xPnicUpgene <- enricher(gene,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.05, qvalueCutoff = 0.05)

ouf <- paste('Pnic_up_KEGG_Enricher.csv',sep ="\t")
# result
write.csv(xPnicUpgene,ouf)
# barplot
barplot(xPnicUpgene,showCategory=40, includeAll=TRUE, cex.names=0.5, font.size = 6)
# bubble plot
dotplot(xPnicUpgene,showCategory=40, font.size = 6)
#enrich plot
xPnicUpgene<-enrichplot::pairwise_termsim(xPnicUpgene)
emapplot(xPnicUpgene,  showCategory = 30, cex_label_category=0.4)
g<-cnetplot(xPnicUpgene, showCategory = 20,cex_label_category = 1, cex_label_gene=0.8)
g + scale_color_manual(values = c("red","blue"))
g1<- cnetplot(xPnicUpgene, node_label="all", showCategory = 30, cex_label_category = 0.8, cex_label_gene=0.45) 
g1 + scale_color_manual(values = c("red","blue"))

p1 <- heatplot(xPnicUpgene)
upset(xPnicUpgene)

options(ggrepel.max.overlaps = Inf) #change to infinate overlap
browseKEGG(xPnicUpgene, "map01100") #Metabolic pathways
browseKEGG(xPnicUpgene, "map01010")
browseKEGG(xPnicUpgene, "ko02010")
browseKEGG(xPnicUpgene, "map00920")
browseKEGG(xPnicUpgene, "ko00920")
browseKEGG(xPnicUpgene, "05415")


library("pathview")
data(xPnicUpgene)
hsa04110 <- pathview(gene.data  = xPnicUpgene,
                     species = "ath",
                     limit  = list(gene=max(abs(geneList)), cpd=1))


