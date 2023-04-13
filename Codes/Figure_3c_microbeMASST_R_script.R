############################################
############################################
####SCRIPT FOR HEATMAP FIGURE 3c############
############################################
############################################


#Setting the working directory and loading packages required 
setwd("D:/POST_DOC/Bile_Acid/microbeMASST_bile_acids/Refined")
library(colorspace)
library(vegan)
library(gplots)
library(ggplot2)

#Loading the input file (available in the source data)
df<-read.csv("all_refined_deltas_heatmap_nogly_notau.csv", check.names = FALSE)

#filtering for bacteria at taxonomic level of phylum and class
df_filter_phylum <- df[,c(2:106, 110)]
df_filter_class <- df[,c(2:106, 109)]

#arranging the row names in ascending order
rnames_phylum <- df_filter_phylum[,106]
mat_phylum <- data.matrix(df_filter_phylum[,1:105])
rownames(mat_phylum) <- rnames_phylum
mat_phylum_2 <- mat_phylum[order(row.names(mat_phylum)), ]
tmat_phylum_2<-t(mat_phylum_2)
write.csv(tmat_phylum_2,"delta_phylum_microbeMASST_filtered.csv")

rnames_class <- df_filter_class[,106]
mat_class <- data.matrix(df_filter_class[,1:105])
rownames(mat_class) <- rnames_class
mat_class_2 <- mat_class[order(row.names(mat_class)), ]
tmat_class_2<-t(mat_class_2)
write.csv(tmat_class_2,"delta_class_microbeMASST_filtered.csv")

#manually remove the delta masses that were observed in all bacterial cultures: 128.96, 128.99, 131.07, 89, 123 and read the corresponding csv for heatmap representation
#manually collapsing all same phylum and class in excel to get presence/absence and reading again 

df_class_collapsed<-read.csv("delta_class_microbeMASST_filtered_collapsed.csv", check.names = FALSE)
rnames_class_collapsed <- df_class_collapsed[,1]
mat_class_collapsed <- data.matrix(df_class_collapsed[,2:18]) #did not keep NA
rownames(mat_class_collapsed) <- rnames_class_collapsed


#heatmap with pheatmap
install.packages("pheatmap")
library(pheatmap)

svg("microbeMASST_class_heatmap.svg")
pheatmap(mat_class_collapsed,cluster_cols = F, cluster_rows = F, main = "microbeMASST_class_level", fontsize_row = 8, fontsize_col = 8, gaps_col = c(2,4,8,9,12,15,16), cellwidth = 15, cellheight = 6, color=colorRampPalette(c("pink", "navy"))(2))
dev.off()

#to make the bar plots 
dfbar<-read.csv("delta_class_microbeMASST_filtered_collapsed_counts_barplot.csv", check.names = FALSE)
library(ggplot2)

dfbar$Class <- factor(df$Class, levels = df$Class) #to avoid automatic reordering
svg("microbeMASST_class_heatmap_barplot.svg")
p<-ggplot(data=dfbar, aes(x=Class, y=Count)) + geom_bar(stat="identity") + theme_classic()
dev.off()


