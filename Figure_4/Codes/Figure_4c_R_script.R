setwd("D:/POST_DOC/Bile_Acid/Reanalysis_public_datasets/MSV000086131_Zoo_animals")

#reading the input file from manual integration in skyline
df<- read.csv("Feature_table_from_skyline_diet.csv", sep =",", quote = "\"")

#reading the file that is manually edited to combine all polyamines into one column
df_2<- read.csv("Feature_table_all_combined_skyline.csv", sep =",", quote = "\"")


####using ggboxplot

library(ggpubr)
library(dplyr)

####all_polyamines_combined
kruskal.test(ATTRIBUTE_diet ~ Polyamine_Area, data = df_2)
pairwise.wilcox.test(df_2$Polyamine_Area, df_2$ATTRIBUTE_diet, p.adjust.method = "BH")
comp <- list( c("Herbivore", "Omnivore"), c("Herbivore", "Carnivore"), c("Omnivore", "Carnivore") )
order <- c("Carnivore", "Omnivore", "Herbivore")
svg("skyline/All_polyamines_diet_stats.svg")
df_2 %>% ggboxplot(x = "ATTRIBUTE_diet", y = "Polyamine_Area", color = "ATTRIBUTE_diet", palette = c("#ECA869", "#B08BBB", "#B5D5C5"), xlab = "All_polyamines", ylab = "Peak area", add = "jitter", add.params = list(color = "ATTRIBUTE_diet", size = 3), order = order) +
  stat_compare_means(comparisons = comp, method = "wilcox")
dev.off()

#zoomed in
svg("skyline/All_polyamines_diet_stats_zoom.svg")
df_2 %>% ggboxplot(x = "ATTRIBUTE_diet", y = "Polyamine_Area", color = "ATTRIBUTE_diet", palette = c("#ECA869", "#B08BBB", "#B5D5C5"), xlab = "All_polyamines", ylab = "Peak area", add = "jitter", add.params = list(color = "ATTRIBUTE_diet", size = 3), order = order, ylim = c(0,20000000)) +
  stat_compare_means(comparisons = comp, method = "wilcox")
dev.off()


#####individual polyamine conjugates - ATTRIBUTE_diet######

####Tri-ornithine
kruskal.test(ATTRIBUTE_diet ~ Tri_orn_area, data = df)
pairwise.wilcox.test(df$Tri_orn_area, df$ATTRIBUTE_diet, p.adjust.method = "BH")
comp <- list( c("Herbivore", "Omnivore"), c("Herbivore", "Carnivore"), c("Omnivore", "Carnivore") )
order <- c("Carnivore", "Omnivore", "Herbivore")
svg("skyline/Tri_ornithine_polyamines_diet_stats.svg")
df %>% ggboxplot(x = "ATTRIBUTE_diet", y = "Tri_orn_area",  color = "ATTRIBUTE_diet", palette = c("#ECA869", "#B08BBB", "#B5D5C5"), xlab = "Trihydroxy_ornithine", ylab = "Peak area",add = "jitter", add.params = list(color = "ATTRIBUTE_diet", size = 3.5), order = order) +
  stat_compare_means(comparisons = comp, method = "wilcox")
dev.off()

####Di-ornithine
kruskal.test(ATTRIBUTE_diet ~ Di_orn_area, data = df)
pairwise.wilcox.test(df$Di_orn_area, df$ATTRIBUTE_diet, p.adjust.method = "BH")
comp <- list( c("Herbivore", "Omnivore"), c("Herbivore", "Carnivore"), c("Omnivore", "Carnivore") )
svg("skyline/Di_ornithine_polyamines_diet_stats.svg")
df %>% ggboxplot(x = "ATTRIBUTE_diet", y = "Di_orn_area", color = "ATTRIBUTE_diet", palette = c("#ECA869", "#B08BBB", "#B5D5C5"), xlab = "Dihydroxy_ornithine", ylab = "Peak area", add = "jitter", add.params = list(color = "ATTRIBUTE_diet", size = 3), order = order) +
  stat_compare_means(comparisons = comp, method = "wilcox")
dev.off()

####Di-acetyl-ornithine
kruskal.test(ATTRIBUTE_diet ~ Di_acetyl_orn_area, data = df)
pairwise.wilcox.test(df_filter$Di_acetyl_orn_area, df_filter$ATTRIBUTE_diet, p.adjust.method = "BH")
comp <- list( c("Herbivore", "Omnivore"), c("Herbivore", "Carnivore"), c("Omnivore", "Carnivore") )
svg("skyline/Di_acetyl_ornithine_polyamines_diet_stats.svg")
df %>% ggboxplot(x = "ATTRIBUTE_diet", y = "Di_acetyl_orn_area", color = "ATTRIBUTE_diet", palette = c("#ECA869", "#B08BBB", "#B5D5C5"), xlab = "Dihydroxy_acetyl_ornithine", ylab = "Peak area", add = "jitter", add.params = list(color = "ATTRIBUTE_diet", size = 3), order = order) +
  stat_compare_means(comparisons = comp, method = "wilcox")
dev.off()
#zoomed in
svg("skyline/Di_acetyl_ornithine_polyamines_diet_stats_zoomed.svg")
df %>% ggboxplot(x = "ATTRIBUTE_diet", y = "Di_acetyl_orn_area", color = "ATTRIBUTE_diet", palette = c("#ECA869", "#B08BBB", "#B5D5C5"), xlab = "Dihydroxy_acetyl_ornithine", ylab = "Peak area", add = "jitter", add.params = list(color = "ATTRIBUTE_diet", size = 3), order = order, ylim = c(0,1000000)) +
  stat_compare_means(comparisons = comp, method = "wilcox")
dev.off()

####Di-putrescine
kruskal.test(ATTRIBUTE_diet ~ Di_put_area, data = df)
pairwise.wilcox.test(df_filter$Di_put_area, df_filter$ATTRIBUTE_diet, p.adjust.method = "BH")
comp <- list( c("Herbivore", "Omnivore"), c("Herbivore", "Carnivore"), c("Omnivore", "Carnivore") )
order <- c("Carnivore", "Omnivore", "Herbivore")
svg("skyline/Di_putrescine_polyamines_diet_stats.svg")
df %>% ggboxplot(x = "ATTRIBUTE_diet", y = "Di_put_area",  color = "ATTRIBUTE_diet", palette = c("#ECA869", "#B08BBB", "#B5D5C5"), xlab = "Dihydroxy_putrescine", ylab = "Peak area", outlier.shape = NA, add = "jitter", add.params = list(color = "ATTRIBUTE_diet", size = 3), order = order) +
  stat_compare_means(comparisons = comp, method = "wilcox")
dev.off()

####Di-cadaverine
kruskal.test(ATTRIBUTE_diet ~ Di_cad_area, data = df)
pairwise.wilcox.test(df_filter$Di_cad_area, df_filter$ATTRIBUTE_diet, p.adjust.method = "BH")
comp <- list( c("Herbivore", "Omnivore"), c("Herbivore", "Carnivore"), c("Omnivore", "Carnivore") )
svg("skyline/Di_cadaverine_polyamines_diet_stats.svg")
df %>% ggboxplot(x = "ATTRIBUTE_diet", y = "Di_cad_area", color = "ATTRIBUTE_diet", palette = c("#ECA869", "#B08BBB", "#B5D5C5"),  xlab = "Dihydroxy_cadaverine", ylab = "Peak area",add = "jitter", add.params = list(color = "ATTRIBUTE_diet", size = 3), order = order) +
  stat_compare_means(comparisons = comp, method = "wilcox")
dev.off()

####Di-acetyl-cadaverine
kruskal.test(ATTRIBUTE_diet ~ Di_acetyl_cad_area, data = df)
pairwise.wilcox.test(df_filter$Di_acetyl_cad_area, df_filter$ATTRIBUTE_diet, p.adjust.method = "BH")
comp <- list( c("Herbivore", "Omnivore"), c("Herbivore", "Carnivore"), c("Omnivore", "Carnivore") )
svg("skyline/Di_acetyl_cadaverine_polyamines_diet_stats.svg")
df %>% ggboxplot(x = "ATTRIBUTE_diet", y = "Di_acetyl_cad_area", color = "ATTRIBUTE_diet", palette = c("#ECA869", "#B08BBB", "#B5D5C5"), xlab = "Dihydroxy_acetyl_cadaverine", ylab = "Peak area", add = "jitter", add.params = list(color = "ATTRIBUTE_diet", size = 3), order = order) +
  stat_compare_means(comparisons = comp, method = "wilcox")
dev.off()

#Di-taurine
kruskal.test(ATTRIBUTE_diet ~ Di_tau_area, data = df)
pairwise.wilcox.test(df$Di_tau_area, df$ATTRIBUTE_diet, p.adjust.method = "BH")
comp <- list( c("Herbivore", "Omnivore"), c("Herbivore", "Carnivore"), c("Omnivore", "Carnivore") )
svg("skyline/11082023_Di_taurine_polyamines_diet_stats.svg")
df %>% ggboxplot(x = "ATTRIBUTE_diet", y = "Di_tau_area", color = "ATTRIBUTE_diet", palette = c("#ECA869", "#B08BBB", "#B5D5C5"), xlab = "Dihydroxy_taurine", ylab = "Peak area", outlier.shape = NA, add = "jitter", add.params = list(color = "ATTRIBUTE_diet", size = 3), order = order) +
  stat_compare_means(comparisons = comp, method = "wilcox")
dev.off()

#Di-glycine
kruskal.test(ATTRIBUTE_diet ~ Di_gly_area, data = df)
pairwise.wilcox.test(df$Di_gly_area, df$ATTRIBUTE_diet, p.adjust.method = "BH")
comp <- list( c("Herbivore", "Omnivore"), c("Herbivore", "Carnivore"), c("Omnivore", "Carnivore") )
svg("skyline/11082023_Di_glycine_polyamines_diet_stats.svg")
df %>% ggboxplot(x = "ATTRIBUTE_diet", y = "Di_gly_area", color = "ATTRIBUTE_diet", palette = c("#ECA869", "#B08BBB", "#B5D5C5"), xlab = "Dihydroxy_glycine", ylab = "Peak area", outlier.shape = NA, add = "jitter", add.params = list(color = "ATTRIBUTE_diet", size = 3), order = order) +
  stat_compare_means(comparisons = comp, method = "wilcox")
dev.off()


