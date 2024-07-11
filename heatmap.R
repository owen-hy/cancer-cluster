## Creating a heatmap
install.packages("pheatmap")
library(pheatmap)

lung_meta40_eda$cluster_id <- k_means_lung$cluster
file_plot <- lung_meta40_eda |>
  arrange(cluster_id)

clust <- file_plot |> 
  select(cluster_id) |>
  mutate(cluster_id = factor(cluster_id))

rownames(clust) = rownames(file_plot)

heat_df <- file_plot[, sapply(file_plot, is.numeric)] |>
  scale() |>
  as.data.frame() |>
  select(-X) |>
  as.matrix()

rownames(heat_df) = rownames(clust)
rownames(clust)

pheatmap(lung_meta40_eda, annotation_row = clust, cluster_cols = FALSE, cluster_rows = FALSE)

??pheatmap

View(file_plot)
View(clust)
View(lung_meta40_eda)