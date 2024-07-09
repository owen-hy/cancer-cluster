hclust(new_lung_meta40)
View(new_lung_meta40)
install.packages("clValid")
library(clValid)

# creating distance matrix (method = euclidean) 
distance_mat <- dist(new_lung_meta40)
# heirarchical clustering 
heir_clust <- hclust(distance_mat, method = 'average')
heir_clust2 <- hclust(distance_mat, method = 'ward.D')
plot(heir_clust)
plot(heir_clust2)
#splitting into clusters 
rect.hclust(heir_clust, k = 3)
rect.hclust(heir_clust2, k = 3)

# K-means clustering
k_means_lung <- eclust(new_lung_meta40, "kmeans", k = 3, nstart = 100)
# before I did it this way 
# k_mean_vals_40 = kmeans(new_lung_meta40, centers = 3, nstart = 100)

# you want the sil_widths to be close to 1, indicates good clustering 
sil <- silhouette(k_means_lung$cluster, dist(new_lung_meta40))
plot(sil)
fviz_silhouette(k_means_lung)

# enhanced heirarchical clustering 
heir_clust_lung <- eclust(new_lung_meta40, "hclust", k = 3, method = "ward.D", graph = FALSE)

sil2 <- silhouette(heir_clust_lung$cluster, dist(new_lung_meta40))
sil2
plot(sil2)
fviz_silhouette(heir_clust_lung)