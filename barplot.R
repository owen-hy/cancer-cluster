## Barplots 

### If you want side by side, replace position = 'fill' with position = 'dodge

current_data <- lung_meta40_hier

current_data <- current_data |> select_if(~!is.numeric(.)) |>
  mutate(stage_numeric = lung_meta40_hier$stage_numeric,
         #Change based on which cluster
         cluster_id = lung_meta40_hier$cluster_id_2) |>
  select(-image_id, -patient_id, -cause_of_death, -stage_at_diagnosis) 

current_data |>
  select(cluster_id, gender) |>
  ggplot(aes(x = factor(cluster_id), fill = gender)) +
  geom_bar(position = 'fill')

current_data |>
  select(cluster_id, mhcII_status) |>
  ggplot(aes(x = factor(cluster_id), fill = mhcII_status)) +
  geom_bar(position = 'fill')

current_data |>
  select(cluster_id, adjuvant_therapy) |>
  ggplot(aes(x = factor(cluster_id), fill = adjuvant_therapy)) +
  geom_bar(position = 'fill')

current_data |>
  select(cluster_id, stage_numeric) |>
  ggplot(aes(x = factor(cluster_id), fill = factor(stage_numeric))) +
  geom_bar(position = 'fill')




  