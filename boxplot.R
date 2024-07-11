## Boxplot for clusters, analysis

## Change according to what you need
## (lung_meta40_hier, lung_meta40_k)
## Smaller cluster is 1 in heir, and 2 in k
current_data <- lung_meta40_hier
# Boxplot of age
current_data |>
  ggplot(aes(x = factor(cluster_id), y = age_at_diagnosis)) +
  geom_boxplot() +
  labs(x = 'Cluster',
       y = 'Age',
       title = 'Age by Cluster')

# Boxplot of pack per years

current_data|>
  ggplot(aes(x = factor(cluster_id), y = pack_years)) +
  geom_boxplot() +
  labs(x = 'Cluster',
       y = 'Pack per Year',
       title = 'Pack per Year by Cluster')

# Boxplot of recurrence day, also accounting for proportion missing

current_data |>
  ggplot(aes(x = factor(cluster_id), y = time_to_recurrence_days)) +
  geom_boxplot() +
  labs(x = 'Cluster',
       y = 'Days till Recurrence',
       title = 'Days till Recurrence by Cluster')

current_data |>
  group_by(cluster_id) |>
  summarize(missing = sum(is.na(time_to_recurrence_days)) / n())


# Boxplot of proportions

current_data |>
  ggplot(aes(x = factor(cluster_id), y = p_CK)) +
  geom_boxplot() +
  labs(x = 'Cluster',
       y = 'Proportion of CK',
       title = 'Proportion of CK by Cluster')

current_data |>
  ggplot(aes(x = factor(cluster_id), y = p_CD8)) +
  geom_boxplot() +
  labs(x = 'Cluster',
       y = 'Proportion of CD8',
       title = 'Proportion of CD8 by Cluster')

current_data |>
  ggplot(aes(x = factor(cluster_id), y = p_CD14)) +
  geom_boxplot() +
  labs(x = 'Cluster',
       y = 'Proportion of CD14',
       title = 'Proportion of CD14 by Cluster')

current_data |>
  ggplot(aes(x = factor(cluster_id), y = p_Other)) +
  geom_boxplot() +
  labs(x = 'Cluster',
       y = 'Proportion of Other Cells',
       title = 'Proportion of Other Cells by Cluster')

current_data |>
  ggplot(aes(x = factor(cluster_id), y = p_CD19)) +
  geom_boxplot() +
  labs(x = 'Cluster',
       y = 'Proportion of CD19',
       title = 'Proportion of CD19 by Cluster')

current_data |>
  ggplot(aes(x = factor(cluster_id), y = p_CD4)) +
  geom_boxplot() +
  labs(x = 'Cluster',
       y = 'Proportion of CD4',
       title = 'Proportion of CD4 by Cluster')
