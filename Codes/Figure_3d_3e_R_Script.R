############################################
############################################
####SCRIPT FOR HEATMAP FIGURE 3d,3e############
############################################
############################################

######Figure 3d#######

#Setting the working directory and loading packages required 
setwd("D:/POST_DOC/Bile_Acid/microbeMASST_bile_acids/Delta_match_in_datasets/Karin_dataset")
df1<-read.csv("Figure_3d_Refined_microbeMASST_Deltas_dihydroxy.csv", check.names = FALSE)
df2<-read.csv("Figure_3d_Refined_microbeMASST_Deltas_tetrahydroxy.csv", check.names = FALSE)
df3<-read.csv("Figure_3d_Refined_microbeMASST_Deltas_trihydroxy.csv", check.names = FALSE)
df4<-read.csv("Figure_3d_Refined_bile_acid_match_Karin_dataset.csv", check.names = FALSE)

#spliting the bile scan column into two
library(tidyr)
df4<-separate(data = df4, col = Bile_acid, into = c("Bile_acid", "Scan"), sep = "_")

#subset for only trihydroxy in karin deltamass
df4_tri<-subset(df4, Bile_acid == 'Trihydroxy')
df4_di<-subset(df4, Bile_acid == 'Dihydroxy')
df4_tetra<-subset(df4, Bile_acid == 'Tetrahydroxy')

df5_tri<-merge(df3,df4_tri,by="delta_mass_round")
df5_di<-merge(df1,df4_di,by="delta_mass_round")
df5_tetra<-merge(df2,df4_tetra,by="delta_mass_round")


#remove duplicated rows
library(dplyr)
df5_tri_unique <- df5_tri %>% distinct()
df5_di_unique <- df5_di %>% distinct()
df5_tetra_unique <- df5_tetra %>% distinct()

df5_tri_unique <- df5_tri[!duplicated(df5_tri[ , c("precursor_mass.x")]), ]
df5_di_unique <- df5_di[!duplicated(df5_di[ , c("precursor_mass.x")]), ]
df5_tetra_unique <- df5_tetra[!duplicated(df5_tetra[ , c("precursor_mass.x")]), ]

df5_tri_unique<-df5_tri_unique[,-c(4:8)]
df5_di_unique<-df5_di_unique[,-c(4:8)]
df5_tetra_unique<-df5_tetra_unique[,-c(4:8)]

df_karin_clusterID<-read.csv("D:/POST_DOC/Bile_Acid/Reanalysis_public_datasets/MSV000080918_Cadaverine/Figure_3d_MSV000080918_cluster_IDs.tsv", sep = "\t", check.names = FALSE)
df_karin_features<-read.csv("D:/POST_DOC/Bile_Acid/Reanalysis_public_datasets/MSV000080918_Cadaverine/Figure_3d_Feature_table.csv", check.names = FALSE)

#Making all column names similar to merge data frame
colnames(df5_tri_unique)[3] <- "precursor_mass"
colnames(df5_di_unique)[3] <- "precursor_mass"
colnames(df5_tetra_unique)[3] <- "precursor_mass"

colnames(df_karin_clusterID)[23] <- "Scan"
colnames(df_karin_clusterID)[28] <- "precursor_mass"

df6_tri<-merge(df5_tri_unique,df_karin_clusterID,by="precursor_mass")
df6_di<-merge(df5_di_unique,df_karin_clusterID,by="precursor_mass")
df6_tetra<-merge(df5_tetra_unique,df_karin_clusterID,by="precursor_mass")

df7_tri<-merge(df6_tri,df_karin_features,by="Scan")
df7_di<-merge(df6_di,df_karin_features,by="Scan")
df7_tetra<-merge(df6_tetra,df_karin_features,by="Scan")

df7_tri<-df7_tri[,-c(1,2,4:32)]
df7_di<-df7_di[,-c(1,2,4:32)]
df7_tetra<-df7_tetra[,-c(1,2,4:32)]

#group by delta mass and aggreagte by sum
df8_tri <- df7_tri %>% 
  group_by(delta_mass_round) %>% 
  summarise(across(everything(), sum))

df8_di <- df7_di %>% 
  group_by(delta_mass_round) %>% 
  summarise(across(everything(), sum))

df8_tetra <- df7_tetra %>% 
  group_by(delta_mass_round) %>% 
  summarise(across(everything(), sum))

write.csv(df8_tri,"microbeMASST_delta_in_karin_dataset_tri_metadata.csv", row.names = FALSE)
write.csv(df8_di,"microbeMASST_delta_in_karin_dataset_di_metadata.csv", row.names = FALSE)
write.csv(df8_tetra,"microbeMASST_delta_in_karin_dataset_tetra_metadata.csv", row.names = FALSE)

#Antibiotic_box_plots
df_box_karin<-read.csv("Antibiotics/Files_for_paper/Figure_3d_microbeMASST_delta_in_karin_dataset_boxplot_tri_antibiotic_final.csv", check.names = FALSE)
df_box_code<-read.csv("Antibiotics/Files_for_paper/Figure_3d_microbeMASST_delta_in_karin_dataset_boxplot_code_tri_final.csv", check.names = FALSE)
library(ggplot2)
svg("microbeMASST_karin_delta_refined_antibiotic_boxplot_tri_X.svg")
ggplot(df_box_karin, aes(x=Antibiotic, y=X, fill=Antibiotic)) +
  geom_jitter(show.legend=FALSE, width=0.25, shape=21) +
  stat_summary(fun = median, show.legend=FALSE, geom="crossbar", color="black", width=0.5, size=0.5) + geom_boxplot(outlier.shape = NA) + theme_classic() + theme(aspect.ratio = 1) 
dev.off()

df_box_karin_di<-read.csv("Antibiotics/Files_for_paper/Figure_3d_microbeMASST_delta_in_karin_dataset_boxplot_di_antibiotic_final.csv", check.names = FALSE)
df_box_code_di<-read.csv("Antibiotics/Files_for_paper/Figure_3d_microbeMASST_delta_in_karin_dataset_boxplot_code_di_final.csv", check.names = FALSE)
library(ggplot2)
svg("microbeMASST_karin_delta_refined_antibiotic_boxplot_di_C.svg")
ggplot(df_box_karin_di, aes(x=Antibiotic, y=C, fill=Antibiotic)) +
  geom_jitter(show.legend=FALSE, width=0.25, shape=21) +
  stat_summary(fun = median, show.legend=FALSE, geom="crossbar", color="black", width=0.5, size=0.5) + geom_boxplot(outlier.shape = NA) + theme_classic() + theme(aspect.ratio = 1) 
dev.off()


df_box_karin_tetra<-read.csv("Antibiotics/Files_for_paper/Figure_3d_microbeMASST_delta_in_karin_dataset_boxplot_tetra_antibiotic_final.csv", check.names = FALSE)
df_box_code_tetra<-read.csv("Antibiotics/Files_for_paper/Figure_3d_microbeMASST_delta_in_karin_dataset_boxplot_code_tetra_final.csv", check.names = FALSE)
library(ggplot2)
svg("microbeMASST_karin_delta_refined_antibiotic_boxplot_tetra_B.svg")
ggplot(df_box_karin_di, aes(x=Antibiotic, y=B, fill=Antibiotic)) +
  geom_jitter(show.legend=FALSE, width=0.25, shape=21) +
  stat_summary(fun = median, show.legend=FALSE, geom="crossbar", color="black", width=0.5, size=0.5) + geom_boxplot(outlier.shape = NA) + theme_classic() + theme(aspect.ratio = 1) 
dev.off()


#Antibiotic_pheatmap with Antibiotic on x-axis and deltas on y -axis

#read the normalized data table (normalized by averaging across yes and no antibiotic)
df_antibiotic<-read.csv("Antibiotics/Figure_3e_Refined_Karin_deltamass_antibiotic_heatmapcsv_all.csv", check.names = FALSE)
rnames <- df_antibiotic[,1]
mat <- data.matrix(df_antibiotic[,2:3])
rownames(mat) <- rnames


#heatmap with pheatmap
library(pheatmap)

svg("microbeMASST_karin_delta_refined_all_antibiotic_final.svg")
pheatmap(mat,cluster_cols = F, cluster_rows = F, main = "microbeMASST_delta_refined_karin", fontsize_row = 12, fontsize_col = 12, cellwidth = 25, cellheight = 20, color=colorRampPalette(c("pink", "navy"))(50))
dev.off()

######Figure 3e#######

#read the normalized data table (normalized by averaging across HFD and NC diet)
df_diet<-read.csv("Diet/Figure_3e_Refined_Karin_deltamass_heatmap_all_diet.csv", check.names = FALSE)
rnames_diet <- df_diet[,1]
mat_diet <- data.matrix(df_diet[,2:3])
rownames(mat_diet) <- rnames_diet


#heatmap with pheatmap
library(pheatmap)

svg("microbeMASST_karin_delta_refined_all_diet_final.svg")
pheatmap(mat_diet,cluster_cols = F, cluster_rows = F, main = "microbeMASST_delta_refined_karin_diet", fontsize_row = 12, fontsize_col = 12, cellwidth = 25, cellheight = 20, color=colorRampPalette(c("pink", "navy"))(50))
dev.off()






