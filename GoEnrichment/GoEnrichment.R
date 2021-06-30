library("clusterProfiler")
# import gene lists
PnicUpgene <- read.csv("/Volumes/WD1/Desktop/laboratory files/Results/Altria project/DEGlist/Wei_reDEGlist/Output/up_Pnic_no rep2.csv",header = T,sep=",")
PnicDowngene<-read.csv("/Volumes/WD1/Desktop/laboratory files/Results/Altria project/DEGlist/Wei_reDEGlist/Output/down_Pnic_no rep2.csv",header = T,sep=",")
FlagUpgene<-read.csv("/Volumes/WD1/Desktop/laboratory files/Results/Altria project/DEGlist/Wei_reDEGlist/Output/up_Flag22_no rep2.csv",header = T,sep=",")
FlagDowngene<-read.csv("/Volumes/WD1/Desktop/laboratory files/Results/Altria project/DEGlist/Wei_reDEGlist/Output/down_Flag22_no rep2.csv",header = T,sep=",")
Flag_PnicUpgene<-read.csv("/Volumes/WD1/Desktop/laboratory files/Results/Altria project/DEGlist/Wei_reDEGlist/Output/up_Flag22_Pnic_no rep2.csv",header = T,sep=",")
Flag_PnicDowngene<-read.csv("/Volumes/WD1/Desktop/laboratory files/Results/Altria project/DEGlist/Wei_reDEGlist/Output/down_Flag22_Pnic_no rep2.csv",header = T,sep=",")
PnicUpgene<-as.factor(PnicUpgene['Row.names'])
PnicDowngene<-as.factor(PnicDowngene['Row.names'])
FlagUpgene<-as.factor(FlagUpgene['Row.names'])
FlagDowngene<-as.factor(FlagDowngene['Row.names'])
Flag_PnicUpgene<-as.factor(Flag_PnicUpgene['Row.names'])
Flag_PnicDowngene<-as.factor(Flag_PnicDowngene['Row.names'])

#import annotation
term2gene <- read.csv("/Users/weiwang/Desktop/GoEnrichment/GOterms annotations.csv",header=T,sep=",")
term2name <- read.csv("/Users/weiwang/Desktop/GoEnrichment/GOterms_GeneList1.csv",header=F,sep=",")

# Enrichment
xPnicUpgene <- enricher(PnicUpgene,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)
ouf <- paste(out_file,sep ="\t")
# result
write.csv(xPnicUpgene,ouf)
# barplot
barplot(x)
# bubble plot
dotplot(x)

xPnicDowngene <- enricher(PnicDowngene,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)
ouf <- paste(out_file,sep ="\t")
# result
write.csv(xPnicDowngene,ouf)
# barplot
barplot(x)
# bubble plot
dotplot(x)

xFlagUpgene <- enricher(FlagUpgene,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)
ouf <- paste(out_file,sep ="\t")
# result
write.csv(xFlagUpgene,ouf)
# barplot
barplot(x)
# bubble plot
dotplot(x)

xFlagDowngene <- enricher(FlagDowngene,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)
ouf <- paste(out_file,sep ="\t")
# result
write.csv(xFlagDowngene,ouf)
# barplot
barplot(x)
# bubble plot
dotplot(x)

xFlag_PnicUpgene <- enricher(Flag_PnicUpgene,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)
ouf <- paste(out_file,sep ="\t")
# result
write.csv(xFlag_PnicUpgene,ouf)
# barplot
barplot(x)
# bubble plot
dotplot(x)

xFlag_PnicDowngene <- enricher(Flag_PnicDowngene,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)
ouf <- paste(out_file,sep ="\t")
# result
write.csv(xFlag_PnicDowngene,ouf)
# barplot
barplot(x)
# bubble plot
dotplot(x)
