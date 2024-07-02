library(tidyverse)
library(MItools)

lungData <- read_rds("./lungData.RDS")
lung_cell <- lungData$lung_cells
lung_meta <- lungData$lung_metadata

unique(lung_cell$pheno) # 6 unqiue cell types (CK, CD8, CD14, Other, CD19, CD4)

# add image-level prevalence of cell types to metadata
lung_meta <- lung_meta |>
  left_join(lung_cell |> 
              group_by(image_id) |>
              summarize(totalCell = n(),
                        p_CK = sum(pheno == "CK")/n(),
                        p_CD8 = sum(pheno == "CD8")/n(),
                        p_CD14 = sum(pheno == "CD14")/n(),
                        p_Other = sum(pheno == "Other")/n(),
                        p_CD19 = sum(pheno == "CD19")/n(),
                        p_CD4 = sum(pheno == "CD4")/n(),),
            join_by(image_id))
                        

# create variables for k function values of each image
lungMeta$k_CK <- NA # k function values for CK
lungMeta$k_CD8 <- NA # k function values for CD8
lungMeta$k_CD14 <- NA # k function values for CD14
lungMeta$k_Other <- NA # k function values for Other
lungMeta$k_CD19 <- NA # k function values for CD19
lungMeta$k_CD4 <- NA # k function values for CD4

# create variables for k cross values of each image
lungMeta$k_CK_CD8 <- NA # k cross values for CK and CD8
lungMeta$k_CK_CD14 <- NA # k cross values for CK and CD14
lungMeta$k_CK_Other <- NA # k cross values for CK and Other
lungMeta$k_CK_CD19 <- NA # k cross values for CK and CD19
lungMeta$k_CK_CD4 <- NA # k cross values for CK and CD4
lungMeta$k_CD8_CD14 <- NA # k cross values for CD8 and CD14
lungMeta$k_CD8_Other <- NA # k cross values for CD8 and Other
lungMeta$k_CD8_CD19 <- NA # k cross values for CD8 and CD19
lungMeta$k_CD8_CD4 <- NA # k cross values for CD8 and CD4
lungMeta$k_CD14_Other <- NA # k cross values for CD14 and Other
lungMeta$k_CD14_CD19 <- NA # k cross values for CD14 and CD19
lungMeta$k_CD14_CD4 <- NA # k cross values for CD14 and CD4
lungMeta$k_Other_CD19 <- NA # k cross values for Other and CD19
lungMeta$k_Other_CD4 <- NA # k cross values for Other and CD4
lungMeta$k_CD19_CD4 <- NA # k cross values for CD19 and CD4