library(pheatmap)
library(readxl)
library(RColorBrewer)
# xls files
options(scipen=999) #prevent scientific notations
options (future.globals.maxSize = 4000 * 1024^7) #increase memory consumption

#vector memory exhausted limit error:
#(from terminal)
#cd ~
#touch .Renviron
#open .Renviron
#add R_MAX_VSIZE=100Gb to  .Renviron file

#Use SSH to access tinkercliff and do the map there
#https://video.vt.edu/media/ARCA+Accessing+clusters+from+the+command+line+via+SSH/1_nkojfb72/176584251
#DOWNLOAD from cload:
#rsync -v wwei6@tinkercliffs2.arc.vt.edu:Pnic_vs_Water.png .

temp<- read_excel("PnicHeatMapData.xlsx", sheet = "Pnic")
heatmap_data <- temp[,-1]
row.names(heatmap_data) <- temp$Gene_ID
summary(c(heatmap_data))

color<-colorRampPalette(c('#436eee','white','#EE0000'))(100)
png(filename = "Pnic_vs_Water.png", height=4000, width=3000, res=500, units="px")
anno_col=data.frame(sampleType=factor(rep(c('Mock','Pnic'), each=21)))
ann_color=list(sampleType=c(mock='#cd0000', normal='#3a5fcd'))
pheatmap(heatmap_data, scale="none", color=color, annotation_col=anno_col, annotation_colors = ann_color,
         fontsize = 8,fontsize_row = 4,show_colnames = F, main = "Pnic vs Mock heatmap",
         legend_breaks=c(4,12),legend_labels = c('low','high'))
dev.off()
