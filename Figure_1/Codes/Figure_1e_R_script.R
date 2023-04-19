############################################
############################################
####SCRIPT FOR HEATMAP FIGURE 1e############
############################################
############################################

#Setting the working directory and loading packages required 
setwd("C:/Users/imoha/OneDrive - University of California, San Diego Health/Bile_Acid/Heatmaps_Fig1")
library(colorspace)
library(vegan)
library(gplots)
library(dplyr)

#Loading the input file (available in the source data)
dfnon<-read.csv("Figure_1e_Non_Refined_counts_morethan2_extracted.csv", check.names = FALSE)

#reshaping the dataframe from long to wide format
dfnon_wide <- reshape(dfnon, idvar = "delta", timevar = "Bile", direction = "wide")
dfnon_wide[is.na(dfnon_wide)] <- 0
dfnon_wide_1<-dfnon_wide[order(dfnon_wide$delta),]
dfnon_wide_1<-dfnon_wide_1[,c(1,4,3,2,7,6,5)]

#rounding the delta masses to 1 to get a bin width of 1Da for the delta masses
dfnon_wide_1$round_delta1<-round(dfnon_wide_1$delta, digits = 1)
df_delta1<-dfnon_wide_1[,-c(9)]

#grouping by deltas and summing over all the counts in non to pentahydroxy columns
df_delta1_grp <- df_delta1 %>%
  group_by(round_delta1) %>%
  summarise(across(everything(), sum))
df_delta1_grp<-df_delta1_grp[,-c(2)]

#input for the blank file of deltas (available in source data) to imput blank values in the delta masses that are not present to ensure a linear scale
df_delta1_blank1<-read.csv("Figure_1e_Delta1_Count_table_Non_Refined_blank.csv", check.names = FALSE)

#subsetting for delta masses that are not in the main delta mass table - basically adding rows with all 0s to make the linear scale
df_delta1_blank1_subset<-subset(df_delta1_blank1, !(round_delta1 %in% df_delta1_grp$round_delta1))

#appending the two dataframes to get the full table
df_delta1_final = rbind(df_delta1_grp,df_delta1_blank1_subset)
df_delta1_final<-as.data.frame(df_delta1_final)

#applying log and storing in a new dataframe 
df_delta1_final_log<-df_delta1_final
df_delta1_final_log[,2:7] <- log10(df_delta1_final_log[2:7] + 1)
df_delta1_final_log2 <- df_delta1_final_log[order(df_delta1_final_log$round_delta1),] #ordering the delta masses in ascending order

#plottinh heatmap using pheatmap package
rnames1 <- df_delta1_final_log2[,1]
matdf1 <- data.matrix(df_delta1_final_log2[,2:7])
rownames(matdf1) <- rnames1
library(pheatmap)
library(viridis)
svg("delta1_nonrefined_heatmap_log.svg")
pheatmap(matdf1,cluster_cols = F, cluster_rows = F, fontsize_row = 2, fontsize_col = 2, cellwidth = 25, cellheight = 0.08, color=viridis(100))
dev.off()

