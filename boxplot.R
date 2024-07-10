## Boxplot for clusters, analysis

## Change according to what you need
current_data <- lung_meta40_eda

# Boxplot of age
current_data |>
  ggplot(aes(x = factor(cluster_id), y = age_at_diagnosis)) +
  geom_boxplot()

# Boxplot of pack per years

current_data|>
  ggplot(aes(x = factor(cluster_id), y = pack_years)) +
  geom_boxplot()

# Boxplot of recurrence day, also accounting for proportion missing

current_data |>
  ggplot(aes(x = factor(cluster_id), y = time_to_recurrence_days)) +
  geom_boxplot()

current_data |>
  group_by(cluster_id) |>
  summarize(missing = sum(is.na(time_to_recurrence_days)) / n())


# Boxplot of proportions

current_data |>
  ggplot(aes(x = factor(cluster_id), y = p_CK)) +
  geom_boxplot()

current_data |>
  ggplot(aes(x = factor(cluster_id), y = p_CD8)) +
  geom_boxplot()

current_data |>
  ggplot(aes(x = factor(cluster_id), y = p_CD14)) +
  geom_boxplot()

current_data |>
  ggplot(aes(x = factor(cluster_id), y = p_Other)) +
  geom_boxplot()

current_data |>
  ggplot(aes(x = factor(cluster_id), y = p_CD19)) +
  geom_boxplot()

current_data |>
  ggplot(aes(x = factor(cluster_id), y = p_CD4)) +
  geom_boxplot()
