library("clusterProfiler")
library(enrichR)
# 导入基因列表

PnicUpgene <- read.csv("/Volumes/WD1/Desktop/laboratory files/Results/Altria project/DEGlist/Wei_reDEGlist/Output/up_Pnic_no rep2.csv",header = T,sep=",")
PnicDowngene<-read.csv("/Volumes/WD1/Desktop/laboratory files/Results/Altria project/DEGlist/Wei_reDEGlist/Output/down_Pnic_no rep2.csv",header = T,sep=",")
FlagUpgene<-read.csv("/Volumes/WD1/Desktop/laboratory files/Results/Altria project/DEGlist/Wei_reDEGlist/Output/up_Flag22_no rep2.csv",header = T,sep=",")
FlagDowngene<-read.csv("/Volumes/WD1/Desktop/laboratory files/Results/Altria project/DEGlist/Wei_reDEGlist/Output/down_Flag22_no rep2.csv",header = T,sep=",")
Flag_PnicUpgene<-read.csv("/Volumes/WD1/Desktop/laboratory files/Results/Altria project/DEGlist/Wei_reDEGlist/Output/up_Flag22_Pnic_no rep2.csv",header = T,sep=",")
Flag_PnicDowngene<-read.csv("/Volumes/WD1/Desktop/laboratory files/Results/Altria project/DEGlist/Wei_reDEGlist/Output/down_Flag22_Pnic_no rep2.csv",header = T,sep=",")
PnicUpgene<-PnicUpgene['Row.names']
PnicDowngene<-PnicDowngene['Row.names']
FlagUpgene<-FlagUpgene['Row.names']
FlagDowngene<-FlagDowngene['Row.names']
Flag_PnicUpgene<-Flag_PnicUpgene['Row.names']
Flag_PnicDowngene<-Flag_PnicDowngene['Row.names']


#import annotation
contents <- read.csv("/Users/weiwang/Desktop/GoEnrichment/Pnic_down_GO.csv",header=T,sep=",")
term2name <- read.csv("/Users/weiwang/Desktop/GoEnrichment/GOterms annotations.csv",header=T,sep=",")
term2gene <- read.csv("/Users/weiwang/Desktop/GoEnrichment/GOterms_GeneList1.csv",header=T,sep=",")
term2gene <- data.frame(term2gene)
term2name <- data.frame(term2name)
head(term2name)
# Enrichment
gene <- as.factor(contents$Gene_ID)
xPnicUpgene <- enricher(gene,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.05, qvalueCutoff = 0.05)

ouf <- paste('Enricher.out.csv',sep ="\t")
# result
write.csv(xPnicUpgene,ouf)
# barplot
barplot(xPnicUpgene)
# bubble plot
dotplot(xPnicUpgene)

xPnicDowngene <- enricher(PnicDowngene,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)
ouf <- paste(out_file,sep ="\t")
# result
write.csv(xPnicDowngene,ouf)
# barplot
barplot(xPnicDowngene)
# bubble plot
dotplot(xPnicDowngene)

xFlagUpgene <- enricher(FlagUpgene,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)
ouf <- paste(out_file,sep ="\t")
# result
write.csv(xFlagUpgene,ouf)
# barplot
barplot(xFlagUpgene)
# bubble plot
dotplot(xFlagUpgene)

xFlagDowngene <- enricher(FlagDowngene,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)
ouf <- paste(out_file,sep ="\t")
# result
write.csv(xFlagDowngene,ouf)
# barplot
barplot(xFlagDowngene)
# bubble plot
dotplot(xFlagDowngene)

xFlag_PnicUpgene <- enricher(Flag_PnicUpgene,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)
ouf <- paste(out_file,sep ="\t")
# result
write.csv(xFlag_PnicUpgene,ouf)
# barplot
barplot(xFlag_PnicUpgene)
# bubble plot
dotplot(xFlag_PnicUpgene)

xFlag_PnicDowngene <- enricher(Flag_PnicDowngene,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)
ouf <- paste(out_file,sep ="\t")
# result
write.csv(xFlag_PnicDowngene,ouf)
# barplot
barplot(xFlag_PnicDowngene)
# bubble plot
dotplot(xFlag_PnicDowngene)
