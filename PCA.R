library(tidyverse)

# Setting up dataset

current_file <- read_csv("./lung_meta40.csv") 

num_file <- current_file |>
  select(image_id, starts_with("p_"), starts_with("k")) |>
  na.omit() 

## Change if using ovarian
image_name <- num_file$image_id

num_mat <- num_file |>
  select(-image_id) |>
  scale() |>
  as.matrix()

rownames(num_mat) <- image_name

# Performing PCA
  
pca <- prcomp(num_mat, scale = T)

pca$rotation

pca$x

# Computing ideal num of PC

pc_var = pca$sdev^2

pve = pc_var / sum(pc_var)
pve

plot(pve, xlab = "Principal Component", ylab = "Proportion of Variance Explained", 
     ylim = c(0, 1), type = "b")
