library(pheatmap)
library(readxl)
library(RColorBrewer)
# xls files
#options(scipen=999) #prevent scientific notations
#options (future.globals.maxSize = 4000 * 1024^7) #increase memory consumption
options(bitmapType='cairo')
#vector memory exhausted limit error:
#(from terminal)
#cd ~
#touch .Renviron
#open .Renviron
#add R_MAX_VSIZE=100Gb to  .Renviron file

#Use SSH to access tinkercliff and do the map there
#https://video.vt.edu/media/ARCA+Accessing+clusters+from+the+command+line+via+SSH/1_nkojfb72/176584251
#DOWNLOAD from ARC:
#scp wwei6@tinkercliffs2.arc.vt.edu:/home/wwei6/Pnic_vs_Water1.png Pnic_vs_Water1.png
#run on ARC code needs to be modified:  

temp<- read_excel("PnicHeatMapData.xlsx", sheet = "Pnic")
temp <- as.data.frame(temp)
rownames(temp) <- temp[, 1]
temp <- temp[, -1]
color<-colorRampPalette(c('#436eee','white','#EE0000'))(100)
png(filename = "Pnic_vs_Water2.png", height=4000, width=3000, res=500, units="px")
anno_col=data.frame(sampleType=factor(c(rep('Mock',21), rep('Pnic',21))))
ann_color=list(sampleType=c(Mock='#cd0000', Pnic='#3a5fcd'))
pheatmap(temp, scale="none", color=color, annotation_col=anno_col, annotation_colors = ann_color,
         fontsize = 4,fontsize_row = 4,show_colnames = F, main = "Pnic vs Mock heatmap",
         legend_breaks=c(4,12),legend_labels = c('low','high'))
dev.off()

#temp<- read_excel("PnicHeatMapData.xlsx", sheet = "Pnic")
#heatmap_data <- temp[,-1]
#row.names(heatmap_data) <- temp$Gene_ID
#summary(c(heatmap_data))

#color<-colorRampPalette(c('#436eee','white','#EE0000'))(100)
#png(filename = "Pnic_vs_Water2.png", height=4000, width=3000, res=500, units="px")
#anno_col=data.frame(sampleType=factor(c(rep('Mock',21), rep('Pnic',21))))
#ann_color=list(sampleType=c(Mock='#cd0000', Pnic='#3a5fcd'))
#pheatmap(heatmap_data, scale="none", color=color, annotation_col=anno_col, annotation_colors = ann_color,
         #fontsize = 4,fontsize_row = 4,show_colnames = F, main = "Pnic vs Mock heatmap",
         #legend_breaks=c(4,12),legend_labels = c('low','high'))
#dev.off()
