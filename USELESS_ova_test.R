
new_ovarian_meta40 = subset(ovarian_meta40, select = -c(sample_id, tma, diagnosis, primary, recurrent, treatment_effect, stage, grade, survival_time, death, BRCA_mutation, age_at_diagnosis, time_to_recurrence, parpi_inhibitor, debulking, totalCell, p_Other, p_Tumor, p_BCell, p_Macrophage, p_THelper, p_CytotoxicT))

new_ovarian_meta40 = na.omit(new_ovarian_meta40)

new_ovarian_meta40 = scale(new_ovarian_meta40)

View(new_ovarian_meta40)


View(ovarian_meta40)

## Other Directions
# We have also begun looking into spectral clustering for non-linear data. Using this method, we determined that we should use four clusters. Using code from [this tutorial](https://rpubs.com/gargeejagtap/SpectralClustering), we have produced some initial results shown below.
```{r}
# Spectral Clustering
# https://rpubs.com/gargeejagtap/SpectralClustering
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
# eigen(laplacian)

# using eigengap to figure out ideal number of cluster
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
```