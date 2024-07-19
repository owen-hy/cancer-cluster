# load libraries
library(dplyr)
library(ggplot2)
library(factoextra)
library(patchwork)
library(cluster)
library(pheatmap)

# Load lung data with Ripley's K calculations using r = 40 microns
lung_meta40 <- read.csv("./lung_meta40.csv")

# Select columns
new_lung_meta40 <- lung_meta40 |>
  select(starts_with("k"), starts_with("p")) |>
  select(-c(pack_years, patient_id))

# Scale data and remove rows containing NA's
new_lung_meta40 <- as.data.frame(scale(new_lung_meta40))
new_lung_meta40 <- na.omit(new_lung_meta40)

# Perform PCA
pca_lung40 <- prcomp(new_lung_meta40)
## View cumulative proportion of variance
summary(pca_lung40)

# Setting a seed
set.seed(42)

# Perform Hierarchical clustering using Ward distance and first 5 principle components
## Produces same dendrogram, but want different objects to draw clusters on
hier_clust2 <- hcut(new_lung_meta40, k = 2)
## 3 clusters
hier_clust3 <- hcut(new_lung_meta40, k = 3)
## 4 clusters
hier_clust4 <- hcut(new_lung_meta40, k = 4)

# View dendrogram and clusters
## 2 clusters
plot(hier_clust2)
rect.hclust(hier_clust2, k = 2)
## 3 clusters
plot(hier_clust3)
rect.hclust(hier_clust3, k = 3)
## 4 clusters
plot(hier_clust4)
rect.hclust(hier_clust4, k = 4)
## silhouette scores for hierarchical clustering
fviz_silhouette(hier_clust2) # 0.11
fviz_silhouette(hier_clust3) # 0.11
fviz_silhouette(hier_clust4) # 0.10

# Create data frame for EDA for hierarchical clustering
lung_meta40_hier <- lung_meta40 |>
  filter(!is.na(k_CK) & !is.na(k_CD8) & !is.na(k_CD14) & !is.na(k_Other) & !is.na(k_Other) & !is.na(k_CD19)
         & !is.na(k_CD4) & !is.na(k_CK_CD8) & !is.na(k_CK_CD14) & !is.na(k_CK_Other) & !is.na(k_CK_CD19) & 
           !is.na(k_CK_CD4) & !is.na(k_CD8_CD14) & !is.na(k_CD8_Other) & !is.na(k_CD8_CD19) & !is.na(k_CD8_CD4)
         & !is.na(k_CD14_Other) & !is.na(k_CD14_CD19) & !is.na(k_CD14_CD4) & !is.na(k_Other_CD19) & 
           !is.na(k_Other_CD4) & !is.na(k_CD19_CD4))

# Add clusters as column
lung_meta40_hier$cluster_id_2 <- hier_clust2$cluster
lung_meta40_hier$cluster_id_3 <- hier_clust3$cluster
lung_meta40_hier$cluster_id_4 <- hier_clust4$cluster

# Count number of patients in each cluster
lung_meta40_hier |>
  group_by(cluster_id_2) |>
  count() # 22 cluster 1, 114 cluster 2

lung_meta40_hier |>
  group_by(cluster_id_3) |>
  count() # 22 cluster 1, 18 cluster 2, 96 cluster 3

lung_meta40_hier |>
  group_by(cluster_id_4) |>
  count() # 19 clutser 1, 18 cluster 2, 3 cluster 3, 96 cluster 4

# create data frame containing medians for all numeric values in each cluster produced using hierarchical clustering
lung_meta40_hier_2 <- lung_meta40_hier |>
  select(-c(cluster_id_3, cluster_id_4))
data<- lung_meta40_hier_2[, sapply(lung_meta40_hier_2, is.numeric)]
cluster_meds_2 <- data |>
  group_by(cluster_id_2) |>
  summarize(across(everything(), median, na.rm = TRUE))

lung_meta40_hier_3 <- lung_meta40_hier |>
  select(-c(cluster_id_2, cluster_id_4))
data<- lung_meta40_hier_3[, sapply(lung_meta40_hier_3, is.numeric)]
cluster_meds_3 <- data |>
  group_by(as.integer(cluster_id_3)) |>
  summarize(across(everything(), median, na.rm = TRUE))

lung_meta40_hier_4 <- lung_meta40_hier |>
  select(-c(cluster_id_3, cluster_id_2))
data <- lung_meta40_hier_4[, sapply(lung_meta40_hier_4, is.numeric)]
cluster_meds_4 <- data |>
  group_by(as.integer(cluster_id_4)) |>
  summarize(across(everything(), median, na.rm = TRUE))

colnames(lung_meta40_hier) <- c("X", "image_id", "patient_id", "gender", "mhcII_status", "age_at_diagnosis",
                                "stage_at_diagnosis", "stage_numeric", "pack_years", "survival_days",
                                "survival_status", "cause_of_death", "adjuvant_therapy",
                                "time_to_recurrence_days", "recurrence_or_lung_ca_death", "total_cell",
                                "p_Tumor", "p_CytotoxicT", "p_Macrophage", "p_Other", "p_BCell", "p_THelper",
                                "k_Tumor", "k_CytotoxicT", "k_Macrophage", "k_Other", "k_BCell", "k_THelper",
                                "k_Tumor_CytotoxicT", "k_Tumor_Macrophage", "k_Tumor_Other", "k_Tumor_BCell",
                                "k_Tumor_THelper", "k_CytotoxicT_Macrophage", "k_CytotoxicT_Other", 
                                "k_CytoxicT_BCell", "k_CytotoxicT_THelper", "k_Macrophage_Other", 
                                "k_Macrophage_BCell", "k_Macrophage_THelper", "k_Other_BCell",
                                "k_Other_THelper", "k_BCell_THelper", "cluster_id_2", "cluster_id_3",
                                "cluster_id_4")

# Making categorical data numeric for heatmaps
lung_meta40_hier2 <- lung_meta40_hier

lung_meta40_hier2$gender <- factor(lung_meta40_hier2$gender, levels = c("M", "F"), labels = c("1","2"))

lung_meta40_hier2$mhcII_status <- factor(lung_meta40_hier2$mhcII_status, levels = c("low", "high"), labels = c("1","2"))

lung_meta40_hier2$adjuvant_therapy <- factor(lung_meta40_hier2$adjuvant_therapy, levels = c("No", "Yes"), labels = c("1","2"))

lung_meta40_hier2$gender = as.integer(lung_meta40_hier2$gender) 
lung_meta40_hier2$mhcII_status = as.integer(lung_meta40_hier2$mhcII_status)
lung_meta40_hier2$adjuvant_therapy = as.integer(lung_meta40_hier2$adjuvant_therapy) 

# Creating heatmaps
## 2 clusters
# set.seed(42)
file_plot <- lung_meta40_hier2 |>
  select(-c(cluster_id_3, cluster_id_4)) |>
  arrange(cluster_id_2) |>
  mutate(new_cluster_id_2 = factor(cluster_id_2)) |>
  select(-cluster_id_2)
clust <- file_plot |> select(new_cluster_id_2)
rownames(clust) = rownames(file_plot)
heat_df <- file_plot[, sapply(file_plot, is.numeric)] |>
  scale() |>
  as.data.frame() |>
  select(-X) |>
  as.matrix()
rownames(heat_df) = rownames(clust)
rownames(clust)
pheatmap(heat_df, 
         annotation_row = clust, 
         cluster_cols = TRUE, 
         cluster_rows = FALSE,
         clustering_method = "ward.D2")
View(file_plot)


## 3 clusters
file_plot <- lung_meta40_hier2 |>
  select(-c(cluster_id_2, cluster_id_4)) |>
  arrange(cluster_id_3) |>
  mutate(new_cluster_id_3 = factor(cluster_id_3)) |>
  select(-cluster_id_3)
clust <- file_plot |> select(new_cluster_id_3)
rownames(clust) = rownames(file_plot)
heat_df <- file_plot[, sapply(file_plot, is.numeric)] |>
  scale() |>
  as.data.frame() |>
  select(-X) |>
  as.matrix()
rownames(heat_df) = rownames(clust)
rownames(clust)
pheatmap(heat_df, 
         annotation_row = clust, 
         cluster_cols = TRUE, 
         cluster_rows = FALSE,
         clustering_method = "ward.D2")
View(file_plot)

## 4 clusters
file_plot <- lung_meta40_hier2 |>
  select(starts_with("k"), starts_with("p"), -pack_years, cluster_id_4) |>
  arrange(cluster_id_4) |>
  mutate(cluster_id = factor(cluster_id_4)) |>
  select(-cluster_id_4)
clust <- file_plot |> select(cluster_id)
rownames(clust) = rownames(file_plot)
heat_df <- file_plot[, sapply(file_plot, is.numeric)] |>
  scale() |>
  as.data.frame() |>
  as.matrix()
rownames(heat_df) = rownames(clust)
rownames(clust)
pheatmap(heat_df, 
         annotation_row = clust, 
         cluster_cols = TRUE, 
         cluster_rows = FALSE,
         clustering_method = "ward.D2")




