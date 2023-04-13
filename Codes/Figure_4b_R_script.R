############################################
############################################
####SCRIPT FOR HEATMAP FIGURE 4b############
############################################
############################################

setwd("D:/POST_DOC/Bile_Acid/fastMASST/Refined")

#reading the fastmasst merged with ReDU table
df1<- read.csv("fastmasst_all_merged_ReDU_refined_final.tsv", sep ="\t", quote = "\"", strip.white=TRUE)
library(dplyr)

#subsetting for disease ontology
df1_subset<- df1[,c(10,21)]

#renaming to maintain consistency
df1_subset[is.na(df1_subset)] = "not specified"
df1_subset$DOIDCommonName[df1_subset$DOIDCommonName=="not applicable"]<-"not specified"
df1_subset$DOIDCommonName[df1_subset$DOIDCommonName=="N/A"]<-"not specified"

#rounding off to two decimal places
df1_subset <- df1_subset %>% mutate_at(vars(delta_mass), funs(round(., 2)))

#grouping by DOIDCommonName
df2_subset<-df1_subset %>%
  group_by(DOIDCommonName, delta_mass) %>%
  mutate(count = n())

#removing duplicates
df3_subset<-df2_subset %>% distinct()

#subsetting for polyamines
df_arg <- subset(df3_subset, delta_mass=='156.10')
df_orn <- subset(df3_subset, delta_mass=='114.08')
df_acetylorn <- subset(df3_subset, delta_mass=='156.09')
df_put <- subset(df3_subset, delta_mass=='70.09')
##df_carbamoylput <- subset(df3_subset, delta_mass=='113.09') -> remove as this maybe Ile/leu conjugated to bile acid
df_acetylput <- subset(df3_subset, delta_mass=='112.10')
df_acetylsperm <- subset (df3_subset, delta_mass=='169.16')

df_cad <- subset(df3_subset, delta_mass=='84.11')
df_acetylcad <- subset(df3_subset, delta_mass=='126.12')
df_DAP <- subset(df3_subset, delta_mass=='56.07')

#concatenating all polyamines files
all_polyamines <- rbind(df_arg, df_orn, df_acetylorn, df_put, df_acetylput, df_acetylsperm, df_cad, df_acetylcad, df_DAP)
all_polyamines<-all_polyamines %>% distinct()

#Renaming diseases to a more controlled vocabulary
all_polyamines <- all_polyamines %>% mutate(DOIDCommonName = recode(DOIDCommonName, "Chagas disease" = "Chagas disease", "circadian rhythm disorders" = "circadian disorders", "sleep deprivation" = "sleep deprivation", "Crohn's disease" = "CD", "ulcerative colitis" = "UC", "inflammatory bowel disease" = "IBD", "obesity" = "obesity","disease NOS" = "disease NOS", "Alzheimer's disease" = "AD", "COVID-19" = "COVID19","acquired immunodeficiency syndrome" = "AIDS","diabetes mellitus" = "diabetes","human immunodeficiency virus infectious disease" = "HIV","Influenza Virus" = "Influenza","Kawasaki diseas" = "Kawasaki disease","mental depression" = "depression","primary bacterial infectious disease" = "Bacterial infection","no disease reported" = "no disease reported","ML import: not available" = "ML import: not available","not collected" = "not collected","not specified" = "not specified" ))

#changing to a wide format
library(tidyr)
all_polyamines_2<-spread(all_polyamines, key = DOIDCommonName, value = count)

all_polyamines_2[is.na(all_polyamines_2)] <- 0

#sub-setting to keep columns of interest
all_polyamines_3<-all_polyamines_2[,-c(5,6)]

#performing normalization by columns
all_polyamines_3[,-1] <- apply(all_polyamines_3[ , -1], 2, function(x){x/sum(x, na.rm=TRUE)})
colnames(all_polyamines_3)[1] <- "disease"

all_polyamines_3<-as.data.frame(all_polyamines_3)

#converting all NaN to 0
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
all_polyamines_3[is.nan(all_polyamines_3)] <- 0

#rearranging columns
all_polyamines_3<-all_polyamines_3[,c(1,4,7,2,3,6,5)]

#removing rows where rowSum is 0
all_polyamines_3 <- all_polyamines_3[rowSums(all_polyamines_3[2:7])>0,]

#converting numerical values to log scale
all_polyamines_4 <- all_polyamines_3[,-c(1)] + 0.0001
all_polyamines_5 <- log(all_polyamines_4)
all_polyamines_6 <- all_polyamines_5*-1

all_polyamines_6['polyamines'] <- all_polyamines_3['disease']

#rearranging columns again
all_polyamines_6<-all_polyamines_6[,c(7,1,2,3,4,5,6)]

library(gplots)
library(ggplot2)

#converting the df to  matrix
rnames <- all_polyamines_6[,1]
mat_all_polyamines_6 <- data.matrix(all_polyamines_6[,2:7])
rownames(mat_all_polyamines_6) <- rnames

#coloring the heatmap based on pink and navy in colorRampPalette in pheatmap
library(RColorBrewer)
breaksList = seq(0, 9.3, by = 0.001)
library(pheatmap)
#scale for this figure: 0 -> 1, 10 -> 0
svg("fastMASST_heatmap_healthcondition_onlypolyamine_log.svg")
pheatmap(mat_all_polyamines_6,cluster_cols = F, cluster_rows = F, color=colorRampPalette(c("navy", "pink"))(length(breaksList)), breaks = breaksList, cellwidth = 15, cellheight = 15)
dev.off()
