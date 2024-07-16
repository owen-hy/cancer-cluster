## Boxplot for clusters, analysis

## Change according to what you need
## (lung_meta40_hier, lung_meta40_k)
## Smaller cluster is 1 in heir
current_data <- lung_meta40_hier

current_data <- current_data[, sapply(current_data, is.numeric)] |>
  select(-X, -stage_numeric, -survival_status, -recurrence_or_lung_ca_death) |>
  ## Change the cluster id if needed
  pivot_longer(cols = -cluster_id_2,
               names_to = "type",
               values_to = "value") 

current_data |>
  ggplot(aes(x = type, y = value, col = factor(cluster_id_2))) +
  geom_boxplot() +
  geom_violin() +
  coord_flip() #This is for all, but you can't really tell much

## Just Proportions

current_data |>
  filter(startsWith(type, "p_")) |>
  ggplot(aes(x = type, y = value, fill = factor(cluster_id_2))) +
  geom_boxplot() +
  coord_flip()

## Just K standalone (except CK, because it was harder to see)

current_data |>
  filter(type %in% c('k_CD8', 'k_CD4', 'k_CD19', 'k_CD14', 'k_Other')) |>
  ggplot(aes(x = type, y = value, fill = factor(cluster_id_2))) +
  geom_boxplot() +
  coord_flip()

current_data |>
  filter(type == 'k_CK') |>
  ggplot(aes(x = type, y = value, fill = factor(cluster_id_2))) +
  geom_boxplot() 

## K Crosses for CK

current_data |>
  filter(startsWith(type, "k_CK_")) |>
  ggplot(aes(x = type, y = value, fill = factor(cluster_id_2))) +
  geom_boxplot() +
  coord_flip()

## K Crosses for CD8

current_data |>
  filter(startsWith(type, "k_CD8_") | endsWith(type, "_CD8"),
         type != 'p_CD8',
         type != 'k_CD8') |>
  ggplot(aes(x = type, y = value, fill = factor(cluster_id_2))) +
  geom_boxplot() +
  coord_flip()

# K Crosses for CD14

current_data |>
  filter(startsWith(type, "k_CD14_")  | endsWith(type, "_CD14"),
         type != 'p_CD14',
         type != 'k_CD14') |>
  ggplot(aes(x = type, y = value, fill = factor(cluster_id_2))) +
  geom_boxplot() +
  coord_flip()

# K Crosses for Other

current_data |>
  filter(startsWith(type, "k_Other_")  | endsWith(type, "_Other"),
         type != 'p_Other',
         type != 'k_Other') |>
  ggplot(aes(x = type, y = value, fill = factor(cluster_id_2))) +
  geom_boxplot() +
  coord_flip()

# Last K Cross for CD19

current_data |>
  filter(startsWith(type, "k_CD19_")  | endsWith(type, "_CD19"),
         type != 'p_CD19',
         type != 'k_CD19') |>
  ggplot(aes(x = type, y = value, fill = factor(cluster_id_2))) +
  geom_boxplot() +
  coord_flip()

# For metadata

current_data |>
  filter(!startsWith(type, "k_"), !startsWith(type, "p_"), 
         type != "pack_years", type != 'age_at_diagnosis') |>
  ggplot(aes(x = type, y = value, fill = factor(cluster_id_2))) +
  geom_boxplot() +
  coord_flip()

current_data |>
  filter(type == "pack_years"| type == 'age_at_diagnosis') |>
  ggplot(aes(x = type, y = value, fill = factor(cluster_id_2))) +
  geom_boxplot() +
  coord_flip()


