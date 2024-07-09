# #hclust(new_lung_meta40)
# #View(new_lung_meta40)
# install.packages("clValid")
# library(clValid)

# creating distance matrix (method = euclidean) 
distance_mat <- dist(new_lung_meta40)
# hierarchical clustering 
hier_clust <- hclust(distance_mat, method = 'ward.D')
plot(hier_clust)
#splitting into clusters 
rect.hclust(heir_clust, k = 3)

# K-means clustering
k_means_lung <- eclust(new_lung_meta40, "kmeans", k = 4, nstart = 100)
# before I did it this way 
# k_mean_vals_40 = kmeans(new_lung_meta40, centers = 3, nstart = 100)

# you want the sil_widths to be close to 1, indicates good clustering 
sil <- silhouette(k_means_lung$cluster, dist(new_lung_meta40))
plot(sil)
fviz_silhouette(k_means_lung)

# enhanced hierarchical clustering 
hier_clust_lung <- eclust(lung_transform3, FUNcluster = "hclust", k = 2, hc_method = "ward.D", graph = TRUE)

sil2 <- silhouette(hier_clust_lung$cluster, dist(new_lung_meta40))
sil2
plot(sil2)
fviz_silhouette(hier_clust_lung)

View(lung_transform)
View(new_lung_meta40)
??eclust

outlierPCOut(new_lung_meta40)

??fviz_nbclust
