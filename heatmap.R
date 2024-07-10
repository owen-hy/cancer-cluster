## Creating a heatmap
library(pheatmap)

file_plot <- lung_meta40_eda_spectral

clust <- file_plot |> select(cluster_id) |>
  mutate(cluster_id = factor(cluster_id))
rownames(clust) = rownames(file_plot)

heat_df <- file_plot[, sapply(file_plot, is.numeric)] |>
  scale() |>
  as.data.frame() |>
  select(-X) |>
  as.matrix()

pheatmap(heat_df, annotation_row = clust)
  
rownames(heat_df) = rownames(clust)
rownames(clust)


