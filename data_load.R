# loading in datasets and necessary libraries
# MItools is downloaded via devtools::install_github("junsoukchoi/MItools")
library(tidyverse)
library(MItools)

lungData <- read_rds("./lungData.RDS")
lung_cell <- lungData$lung_cells
lung_meta <- lungData$lung_metadata
lung_combined <- merge(lungData$lung_cells, lungData$lung_metadata)

ovarianData <- read_rds("./ovarianData.RDS")
ovarian_cell <- ovarianData$ovarian_cells
ovarian_meta <- ovarianData$ovarian_metadata
ovarian_combined <- merge(ovarianData$ovarian_cells, ovarianData$ovarian_metadata)

lung_filtered <- lung_combined %>% 
  group_by(patient_id, image_id) %>% 
  count() %>%
  group_by(patient_id) %>% 
  arrange(desc(n)) %>% 
  distinct(patient_id, .keep_all = TRUE) %>% 
  left_join(lung_combined[ , c("image_id" , "slide_id", "x", "y", "pheno", "gender", "mhcII_status", "age_at_diagnosis", "stage_numeric", "pack_years", "survival_days", "cause_of_death", "adjuvant_therapy", "time_to_recurrence_days", "recurrence_or_lung_ca_death")], by = "image_id")

unique(lung_filtered$image_id)
## Calculating K and K cross
