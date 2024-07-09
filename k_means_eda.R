library(dplyr)
library(ggplot2)
library(factoextra)
library(patchwork) #package to plot multiple graphs at the same time 
names(lung_meta)

lung_meta40 <- read.csv("lung_meta40.csv")
lung_meta30 <- read.csv("lung_meta30.csv")
lung_meta20 <- read.csv("lung_meta20.csv")

ovarian_meta40 <- read.csv("ovarian_meta40.csv")
ovarian_meta30 <- read.csv("ovarian_meta30.csv")
ovarian_meta20 <- read.csv("ovarian_meta20.csv")

#removing all columns except the k-value ones 
new_lung_meta40 = subset(lung_meta40, select = -c(image_id, patient_id, gender, mhcII_status, age_at_diagnosis, stage_at_diagnosis, stage_numeric, pack_years, survival_days, survival_status, cause_of_death, adjuvant_therapy, time_to_recurrence_days, recurrence_or_lung_ca_death, totalCell))

new_lung_meta30 = subset(lung_meta30, select = -c(image_id, patient_id, gender, mhcII_status, age_at_diagnosis, stage_at_diagnosis, stage_numeric, pack_years, survival_days, survival_status, cause_of_death, adjuvant_therapy, time_to_recurrence_days, recurrence_or_lung_ca_death, totalCell))

new_lung_meta20 = subset(lung_meta20, select = -c(image_id, patient_id, gender, mhcII_status, age_at_diagnosis, stage_at_diagnosis, stage_numeric, pack_years, survival_days, survival_status, cause_of_death, adjuvant_therapy, time_to_recurrence_days, recurrence_or_lung_ca_death, totalCell))

# removing all the na entries 
new_lung_meta40 = na.omit(new_lung_meta40)
new_lung_meta30 = na.omit(new_lung_meta30)
new_lung_meta20 = na.omit(new_lung_meta20)

# scaling the data 
new_lung_meta40 = scale(new_lung_meta40)
new_lung_meta30 = scale(new_lung_meta30)
new_lung_meta20 = scale(new_lung_meta20)

# determining optimal number of clusters
wss_lung = fviz_nbclust(new_lung_meta40, kmeans, method = "wss")
silhouette_lung = fviz_nbclust(new_lung_meta40, kmeans, method = "silhouette")
gap_stat_lung = fviz_nbclust(new_lung_meta40, kmeans, method = "gap_stat")
wss_lung + silhouette_lung + gap_stat_lung

# n = 10 # max number of clusters to consider 
# wss = numeric(n) # the within-cluster sum of squares for each number of clusters   
# for (i in 1:n) {
#   km.out <- kmeans(new_lung_meta40, centers = i, nstart = 100)
#   wss[i] <- km.out$tot.withinss
# }
# wss_df <- tibble(clusters = 1:n, wss = wss)
# scree_plot_lung40 <- ggplot(wss_df, aes(x = clusters, y = wss)) + geom_line() + geom_point() + scale_x_continuous(breaks=c(1:10)) + ggtitle("Elbow plot 40 microns for lung data")
# scree_plot_lung40
# 
# n = 10 
# wss = numeric(n) 
# for (i in 1:n) {
#   km.out <- kmeans(new_lung_meta30, centers = i, nstart = 100)
#   wss[i] <- km.out$tot.withinss
# }
# wss_df <- tibble(clusters = 1:n, wss = wss)
# scree_plot_lung30 <- ggplot(wss_df, aes(x = clusters, y = wss)) + geom_line() + geom_point() + scale_x_continuous(breaks=c(1:10)) + ggtitle("Elbow plot 30 microns for lung data")
# scree_plot_lung30
# 
# 
# n = 10 
# wss = numeric(n) 
# for (i in 1:n) {
#   km.out <- kmeans(new_lung_meta20, centers = i, nstart = 100)
#   wss[i] <- km.out$tot.withinss
# }
# wss_df <- tibble(clusters = 1:n, wss = wss)
# scree_plot_lung20 <- ggplot(wss_df, aes(x = clusters, y = wss)) + geom_line() + geom_point() + scale_x_continuous(breaks=c(1:10)) + ggtitle("Elbow plot 20 microns for lung data")
# scree_plot_lung20
# 
# scree_plot_lung40 + scree_plot_lung30 + scree_plot_lung20

# visualizing clusters 
k_mean_vals_40 = kmeans(new_lung_meta40, centers = 5, nstart = 100)
lung_cluster_40 = fviz_cluster(k_mean_vals, data = new_lung_meta40, 
               geom = "point",
               ellipse.type = "convex", 
               ggtheme = theme_bw(), main = "Cluster plot 40 microns for lung")

k_mean_vals_30 = kmeans(new_lung_meta30, centers = 5, nstart = 100)
lung_cluster_30 = fviz_cluster(k_mean_vals, data = new_lung_meta30, 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw(), main = "Cluster plot 30 microns for lung")

k_mean_vals_20 = kmeans(new_lung_meta20, centers = 5, nstart = 100)
lung_cluster_20 = fviz_cluster(k_mean_vals, data = new_lung_meta20, 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw(), main = "Cluster plot 20 microns for lung")

lung_cluster_40 + lung_cluster_30 + lung_cluster_20

# now doing the same for ovarian data 
names(ovarian_meta40)

new_ovarian_meta40 = subset(ovarian_meta40, select = -c(sample_id, tma, diagnosis, primary, recurrent, treatment_effect, stage, grade, survival_time, death, BRCA_mutation, age_at_diagnosis, time_to_recurrence, parpi_inhibitor, debulking, totalCell, p_Other, p_Tumor, p_BCell, p_Macrophage, p_THelper, p_CytotoxicT))

new_ovarian_meta30 = subset(ovarian_meta30, select = -c(sample_id, tma, diagnosis, primary, recurrent, treatment_effect, stage, grade, survival_time, death, BRCA_mutation, age_at_diagnosis, time_to_recurrence, parpi_inhibitor, debulking, totalCell, p_Other, p_Tumor, p_BCell, p_Macrophage, p_THelper, p_CytotoxicT))

new_ovarian_meta20 = subset(ovarian_meta20, select = -c(sample_id, tma, diagnosis, primary, recurrent, treatment_effect, stage, grade, survival_time, death, BRCA_mutation, age_at_diagnosis, time_to_recurrence, parpi_inhibitor, debulking, totalCell, p_Other, p_Tumor, p_BCell, p_Macrophage, p_THelper, p_CytotoxicT))


new_ovarian_meta40 = na.omit(new_ovarian_meta40)
new_ovarian_meta30 = na.omit(new_ovarian_meta30)
new_ovarian_meta20 = na.omit(new_ovarian_meta20)

new_ovarian_meta40 = scale(new_ovarian_meta40)
new_ovarian_meta30 = scale(new_ovarian_meta30)
new_ovarian_meta20 = scale(new_ovarian_meta20)


# n = 10 
# wss = numeric(n)
# for (i in 1:n) {
#   km.out <- kmeans(new_ovarian_meta40, centers = i, nstart = 100)
#   wss[i] <- km.out$tot.withinss
# }
# wss_df <- tibble(clusters = 1:n, wss = wss)
# scree_plot_ovarian40 <- ggplot(wss_df, aes(x = clusters, y = wss)) + geom_line() + geom_point() + scale_x_continuous(breaks=c(1:10)) + ggtitle("Elbow plot 40 microns for ovarian data")
# scree_plot_ovarian40
# 
# 
# wss = numeric(n)
# for (i in 1:n) {
#   km.out <- kmeans(new_ovarian_meta30, centers = i, nstart = 100)
#   wss[i] <- km.out$tot.withinss
# }
# wss_df <- tibble(clusters = 1:n, wss = wss)
# scree_plot_ovarian30 <- ggplot(wss_df, aes(x = clusters, y = wss)) + geom_line() + geom_point() + scale_x_continuous(breaks=c(1:10)) + ggtitle("Elbow plot 30 microns for ovarian data")
# scree_plot_ovarian30
# 
# 
# wss = numeric(n)
# for (i in 1:n) {
#   km.out <- kmeans(new_ovarian_meta20, centers = i, nstart = 100)
#   wss[i] <- km.out$tot.withinss
# }
# wss_df <- tibble(clusters = 1:n, wss = wss)
# scree_plot_ovarian20 <- ggplot(wss_df, aes(x = clusters, y = wss)) + geom_line() + geom_point() + scale_x_continuous(breaks=c(1:10)) + ggtitle("Elbow plot 20 microns for ovarian data")
# scree_plot_ovarian20
# 
# scree_plot_ovarian40 + scree_plot_ovarian30 + scree_plot_ovarian20


# width sum of squares method (elbow plot)
wss_ova = fviz_nbclust(new_ovarian_meta40, kmeans, method = "wss")

# silhouette method
silhouette_ova = fviz_nbclust(new_ovarian_meta40, kmeans, method = "silhouette")

# gap statistic method 
gap_stat_ova = fviz_nbclust(new_ovarian_meta40, kmeans, method = "gap_stat")

wss_ova + silhouette_ova + gap_stat_ova

