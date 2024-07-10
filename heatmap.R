## Creating a heatmap
install.packages("pheatmap")
library(pheatmap)

file_plot <- lung_meta40_eda_spectral |>
  arrange(cluster_id)

clust <- file_plot |> select(cluster_id) |>
  mutate(cluster_id = factor(cluster_id))
rownames(clust) = rownames(file_plot)

heat_df <- file_plot[, sapply(file_plot, is.numeric)] |>
  scale() |>
  select(-X) |>
  as.matrix()

rownames(heat_df) = rownames(clust)
rownames(clust)

pheatmap(heat_df, annotation_row = clust, cluster_cols = TRUE, cluster_rows = TRUE)
  

View(file_plot)

