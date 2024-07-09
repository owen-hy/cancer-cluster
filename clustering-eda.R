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
  pca_lung40$x[,1:5],
  FUNcluster = hcut,
  method = c("wss"),
  diss = NULL,
  k.max = 10,
  nboot = 100,
  verbose = interactive(),
  barfill = "steelblue",
  barcolor = "steelblue",
  linecolor = "steelblue",
  print.summary = TRUE
)

hier_clust <- hcut(pca_lung40$x[,1:5], k = 2)
plot(hier_clust)
rect.hclust(hier_clust, k = 2)
  
  

fviz_silhouette((hier_clust))

kmeans_lung40 <- kmeans(data.frame(pca_lung40$x[, 1:5]), centers = 2)
kmeans_lung40_e <- eclust(pca_lung40$x[, 1:5], "kmeans", k = 2, nstart = 100)
fviz_pca_ind(pca_lung40, habillage = kmeans_lung40$cluster, addEllipses = TRUE)
fviz_silhouette(kmeans_lung40_e)

lung_meta40_eda <- lung_meta40 |>
  filter(!is.na(p_CK) & !is.na(p_CD8) & !is.na(p_CD14) & !is.na(p_Other) & !is.na(p_CD19) & !is.na(p_CD4)
         & !is.na(k_CK) & !is.na(k_CD8) & !is.na(k_CD14) & !is.na(k_Other) & !is.na(k_Other) & !is.na(k_CD19)
         & !is.na(k_CD4) & !is.na(k_CK_CD8) & !is.na(k_CK_CD14) & !is.na(k_CK_Other) & !is.na(k_CK_CD19) & 
           !is.na(k_CK_CD4) & !is.na(k_CD8_CD14) & !is.na(k_CD8_Other) & !is.na(k_CD8_CD19) & !is.na(k_CD8_CD4)
         & !is.na(k_CD14_Other) & !is.na(k_CD14_CD19) & !is.na(k_CD14_CD4) & !is.na(k_Other_CD19) & 
           !is.na(k_Other_CD4) & !is.na(k_CD19_CD4))

lung_meta40_eda$cluster_id <- hier_clust$cluster
lung_meta40_eda |>
  group_by(cluster_id) |>
  count()

data <- lung_meta40_eda[, sapply(lung_meta40_eda, is.numeric)]
cluster_meds <- data |>
  group_by(cluster_id) |>
  summarize(across(everything(), median, na.rm = TRUE))

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


# Spectral Clustering

## Create Similarity matrix of Euclidean distance between points
S <- as.matrix(dist(pca_lung40$x[,1:5]))


## Create Degree matrix
D <- matrix(0, nrow=nrow(pca_lung40$x[,1:5]), ncol = nrow(pca_lung40$x[,1:5])) # empty nxn matrix

for (i in 1:nrow(pca_lung40$x[,1:5])) {
  
  # Find top 10 nearest neighbors using Euclidean distance
  index <- order(S[i,])[2:11]
  
  # Assign value to neighbors
  D[i,][index] <- 1 
}

# find mutual neighbors
D = D + t(D) 
D[ D == 2 ] = 1

# find degrees of vertices
degrees = colSums(D) 
n = nrow(D)


## Compute Laplacian matrix
# Since k > 2 clusters (3), we normalize the Laplacian matrix:
laplacian = ( diag(n) - diag(degrees^(-1/2)) %*% D %*% diag(degrees^(-1/2)) )
eigen(laplacian)

which.max(diff(sort(eigen(laplacian)$values, decreasing = F))) + 1


## Compute eigenvectors
eigenvectors = eigen(laplacian, symmetric = TRUE)
n = nrow(laplacian)
eigenvectors = eigenvectors$vectors[,(n - 2):(n - 1)]


set.seed(1748)
## Run Kmeans on eigenvectors
sc = kmeans(eigenvectors, 4)


## Pull clustering results
sc_results = cbind(pca_lung40$x[,1:5], cluster = as.factor(sc$cluster))
head(sc_results)

as.data.frame(sc_results) |>
  group_by(cluster) |>
  count()

ggplot(data = sc_results, aes(x = PC1, y = PC2)) + aes(color = factor(cluster)) + geom_point()

lung_meta40_eda_spectral <- lung_meta40 |>
  filter(!is.na(p_CK) & !is.na(p_CD8) & !is.na(p_CD14) & !is.na(p_Other) & !is.na(p_CD19) & !is.na(p_CD4)
         & !is.na(k_CK) & !is.na(k_CD8) & !is.na(k_CD14) & !is.na(k_Other) & !is.na(k_Other) & !is.na(k_CD19)
         & !is.na(k_CD4) & !is.na(k_CK_CD8) & !is.na(k_CK_CD14) & !is.na(k_CK_Other) & !is.na(k_CK_CD19) & 
           !is.na(k_CK_CD4) & !is.na(k_CD8_CD14) & !is.na(k_CD8_Other) & !is.na(k_CD8_CD19) & !is.na(k_CD8_CD4)
         & !is.na(k_CD14_Other) & !is.na(k_CD14_CD19) & !is.na(k_CD14_CD4) & !is.na(k_Other_CD19) & 
           !is.na(k_Other_CD4) & !is.na(k_CD19_CD4))

clusters <- sc_results[, 6]

lung_meta40_eda_spectral$cluster_id <- clusters
lung_meta40_eda_spectral |>
  group_by(cluster_id) |>
  count()

data <- lung_meta40_eda_spectral[, sapply(lung_meta40_eda_spectral, is.numeric)]
cluster_meds_spectral <- data |>
  group_by(cluster_id) |>
  summarize(across(everything(), median, na.rm = TRUE))

fviz_silhouette(silhouette(sc$cluster, dist(eigenvectors)))
