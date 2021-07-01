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

temp<- read_excel("PnicHeatMapData.xlsx", sheet = "GOI")
temp <- data.frame(temp)

rownames(temp) <- temp[, 1]
temp <- temp[, -1]
data_subset <- as.matrix(temp)


#cat_df = data.frame("category" = c(rep('Mock',21), rep('Pnic',21)))


color<-colorRampPalette(c('#436eee','white','#EE0000'))(100)
anno_col=data.frame(sampleType=factor(c(rep('Mock',21), rep('Pnic',21), rep('Flag22', 20), rep('Flag22_Pnic', 21))))
rownames(anno_col)=colnames(data_subset)
ann_color=list(sampleType=c(Mock='#cd0000', Pnic='#3a5fcd', Flag22='#00FF00', Flag22_Pnic='#ffa500'))

All_my_heatmap <- pheatmap(data_subset, scale="none", color=color, annotation_col=anno_col, annotation_colors = ann_color, show_colnames=T, cluster_rows = F, main = "Treatments vs Mock heatmap",legend_breaks=c(4,12),fontsize=8, fontsize_row=6, fontsize_col=4)
save_pheatmap_png <- function(x, filename, width=4800, height=4000, res = 700) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(All_my_heatmap, "all_my_heatmap.png")

#Pnic vs Mock
anno_col=data.frame(sampleType=factor(c(rep('Mock',21), rep('Pnic',21))))
ann_color=list(sampleType=c(Mock='#cd0000', Pnic='#3a5fcd'))

data<-data_subset[c('T1NCMock1':'T1NCpath1','T1NCpath1':'T1NCflag1')]

Pnic_my_heatmap_noCluster <- pheatmap(data_subset[], scale="none", color=color, annotation_col=anno_col, annotation_colors = ann_color, show_colnames=T, cluster_rows = F, main = "Pnic vs Mock heatmap",legend_breaks=c(4,12),fontsize=8, fontsize_row=6, fontsize_col=4)
Pnic_my_heatmap_Cluster<-pheatmap(data_subset, scale="none", color=color, annotation_col=anno_col, annotation_colors = ann_color, show_colnames=T, cluster_rows = T, main = "Pnic vs Mock heatmap",legend_breaks=c(4,12),fontsize=8, fontsize_row=6, fontsize_col=4)


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
