library(dplyr)
library(ggplot2)
library(factoextra)
names(lung_meta)

#removing all columns except the k-value ones 
lung_meta2 = subset(lung_meta, select = -c(image_id, patient_id, gender, mhcII_status, age_at_diagnosis, stage_at_diagnosis, stage_numeric, pack_years, survival_days, survival_status, cause_of_death, adjuvant_therapy, time_to_recurrence_days, recurrence_or_lung_ca_death, totalCell) )
View(lung_meta2)
# removing all the na entries 
lung_meta2 = na.omit(lung_meta2)
# scaling the data 
lung_meta2 = scale(lung_meta2)


n = 10 # max number of clusters to consider 
wss = numeric(n) # the within-cluster sum of squares for each number of clusters   
for (i in 1:n) {
  km.out <- kmeans(lung_meta2, centers = i, nstart = 100)
  wss[i] <- km.out$tot.withinss
}

k_mean_vals = kmeans(lung_meta2, centers = 5, nstart = 100)
fviz_cluster(k_mean_vals, data = lung_meta2, 
               geom = "point",
               ellipse.type = "convex", 
               ggtheme = theme_bw())
