library(factoextra)

lung_transform3 <- as.data.frame(-pca$x[,1:3])

lung_transform <- as.data.frame(-pca$x[,1:8])

fviz_nbclust(lung_transform, FUNcluster = kmeans, method = 'wss')

k_means_pca <- kmeans(lung_transform, centers = 5, nstart = 25)

current_file$cluster <- k_means_pca$cluster

plot_mat <- as.data.frame(cbind(num_mat, k_means_pca$cluster))|>
  mutate(V28 = factor(V28)) |>
  rename('cluster' = 'V28')

autoplot(pca, data = plot_mat, color = 'cluster')




