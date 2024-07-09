library(dplyr)
library(ggplot2)
library(factoextra)
library(cluster)
# library(patchwork) #package to plot multiple graphs at the same time 
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
new_lung_meta40 = as.data.frame(scale(new_lung_meta40)) |>
  select(-c(X, p_CK, p_CD8, p_CD14, p_Other, p_CD19, p_CD4)) 

new_lung_meta30 = as.data.frame(scale(new_lung_meta30)) |>
  select(-X)
new_lung_meta20 = as.data.frame(scale(new_lung_meta20)) |>
  select(-X)

set.seed(42)

# lung metadata with k computed at r = 40 micorns
pca_lung40 <- prcomp(new_lung_meta40)
summary(pca_lung40)

fviz_nbclust(
  pca_lung40$x[,1:7],
  FUNcluster = hcut,
  method = c("silhouette"),
  diss = NULL,
  k.max = 10,
  nboot = 100,
  verbose = interactive(),
  barfill = "steelblue",
  barcolor = "steelblue",
  linecolor = "steelblue",
  print.summary = TRUE
)

kmeans_lung40 <- kmeans(data.frame(pca_lung40$x[, 1:7]), centers = 2)
fviz_pca_ind(pca_lung40, habillage = kmeans_lung40$cluster, addEllipses = TRUE)

lung_meta40_eda <- lung_meta40 |>
  filter(!is.na(p_CK) & !is.na(p_CD8) & !is.na(p_CD14) & !is.na(p_Other) & !is.na(p_CD19) & !is.na(p_CD4)
         & !is.na(k_CK) & !is.na(k_CD8) & !is.na(k_CD14) & !is.na(k_Other) & !is.na(k_Other) & !is.na(k_CD19)
         & !is.na(k_CD4) & !is.na(k_CK_CD8) & !is.na(k_CK_CD14) & !is.na(k_CK_Other) & !is.na(k_CK_CD19) & 
           !is.na(k_CK_CD4) & !is.na(k_CD8_CD14) & !is.na(k_CD8_Other) & !is.na(k_CD8_CD19) & !is.na(k_CD8_CD4)
         & !is.na(k_CD14_Other) & !is.na(k_CD14_CD19) & !is.na(k_CD14_CD4) & !is.na(k_Other_CD19) & 
           !is.na(k_Other_CD4) & !is.na(k_CD19_CD4))

lung_meta40_eda$cluster_id <- kmeans_lung40$cluster

# Summary statistics
## Number of patients in each cluster
lung_meta40_eda |>
  count(cluster_id) # 42 patients in cluster 1, 94 patients in cluster 2

lung_meta40_eda |>
  group_by(cluster_id) |>
  summarize(mean_survival = mean(survival_days, na.rm = TRUE), 
            median_survival = median(survival_days, na.rm = TRUE),
            sd_survival = sd(survival_days, na.rm = TRUE),
            mean_pack = mean(pack_years, na.rm = TRUE), 
            median_pack = median(pack_years, na.rm = TRUE),
            sd_pack = sd(pack_years, na.rm = TRUE),
            mean_age = mean(age_at_diagnosis, na.rm = TRUE), 
            median_age = median(age_at_diagnosis, na.rm = TRUE),
            sd_age = sd(age_at_diagnosis, na.rm = TRUE))

lung_meta40_eda |>
  group_by(cluster_id) |>
  summarize(prop_stage1 = sum(stage_numeric == 1)/n(),
            prop_stage2 = sum(stage_numeric == 2)/n(),
            prop_stage3 = sum(stage_numeric == 3)/n(),
            prop_stage4 = sum(stage_numeric == 4)/n()) # higher proportion of stage 2 patients in cluster 2

