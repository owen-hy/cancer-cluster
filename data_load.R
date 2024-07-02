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

## Calculating K and K cross
