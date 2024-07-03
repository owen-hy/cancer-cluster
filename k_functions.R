library(tidyverse)
library(MItools)

lungData <- read_rds("./lungData.RDS")
lung_cell <- lungData$lung_cells
lung_meta <- lungData$lung_metadata

cell_types <- unique(lung_cell$pheno) # 6 unqiue cell types (CK, CD8, CD14, Other, CD19, CD4)

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
lung_meta$k_CK <- NA # k function values for CK
lung_meta$k_CD8 <- NA # k function values for CD8
lung_meta$k_CD14 <- NA # k function values for CD14
lung_meta$k_Other <- NA # k function values for Other
lung_meta$k_CD19 <- NA # k function values for CD19
lung_meta$k_CD4 <- NA # k function values for CD4

# create variables for k cross values of each image
lung_meta$k_CK_CD8 <- NA # k cross values for CK and CD8
lung_meta$k_CK_CD14 <- NA # k cross values for CK and CD14
lung_meta$k_CK_Other <- NA # k cross values for CK and Other
lung_meta$k_CK_CD19 <- NA # k cross values for CK and CD19
lung_meta$k_CK_CD4 <- NA # k cross values for CK and CD4
lung_meta$k_CD8_CD14 <- NA # k cross values for CD8 and CD14
lung_meta$k_CD8_Other <- NA # k cross values for CD8 and Other
lung_meta$k_CD8_CD19 <- NA # k cross values for CD8 and CD19
lung_meta$k_CD8_CD4 <- NA # k cross values for CD8 and CD4
lung_meta$k_CD14_Other <- NA # k cross values for CD14 and Other
lung_meta$k_CD14_CD19 <- NA # k cross values for CD14 and CD19
lung_meta$k_CD14_CD4 <- NA # k cross values for CD14 and CD4
lung_meta$k_Other_CD19 <- NA # k cross values for Other and CD19
lung_meta$k_Other_CD4 <- NA # k cross values for Other and CD4
lung_meta$k_CD19_CD4 <- NA # k cross values for CD19 and CD4

# for loop to iterate through each image in the metadata
for(i in 1:nrow(lung_meta)){
  # temporary dataset for image "i" in the ith iteration of this "for" loop
  tempData <- lungData$lung_cells |>
    filter(image_id == lung_meta$image_id[i])

  # point pattern object for all of image i
  pppObj <- spatstat.geom::ppp(
    x = tempData$x, y = tempData$y,
    xrange = range(tempData$x), yrange = range(tempData$y),
    marks = factor(tempData$pheno))

  if(lung_meta$p_CK[i] > 0){
    # k function estimated for image i tumor cells using isotropic correction
    kEstCK_temp <- permute.envelope(pppObj, fun = Kest, funargs = list(i = "CK", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of 30 microns
    lung_meta$k_CK[i] <- (kEstCK_temp$obs - kEstCK_temp$mmean)[
      which.min(abs(kEstCK_temp$r - 30))]
  }
  else{
    # if there are no tumor cells, record NA for this image's k function
    lung_meta$k_CK[i] <- NA
  }

  # CD8
  if(lung_meta$p_CD8[i] > 0){
    kEstCD8_temp <- permute.envelope(pppObj, fun = Kest, funargs = list(i = "CD8", correction = "isotropic"))
    lung_meta$k_CD8[i] <- (kEstCD8_temp$obs - kEstCD8_temp$mmean)[
      which.min(abs(kEstCD8_temp$r - 30))]
  }
  else{
    lung_meta$k_CD8[i] <- NA
  }

  # # CD14
  if(lung_meta$p_CD14[i] > 0){
    kEstCD14_temp <- permute.envelope(pppObj, fun = Kest, funargs = list(i = "CD14", correction = "isotropic"))
    lung_meta$k_CD14[i] <- (kEstCD14_temp$obs - kEstCD14_temp$mmean)[
      which.min(abs(kEstCD14_temp$r - 30))]
  }
  else{
    lung_meta$k_CD14[i] <- NA
  }

  # # Other
  if(lung_meta$p_Other[i] > 0){
    kEstOther_temp <- permute.envelope(pppObj, fun = Kest, funargs = list(i = "Other", correction = "isotropic"))
    lung_meta$k_Other[i] <- (kEstOther_temp$obs - kEstOther_temp$mmean)[
      which.min(abs(kEstOther_temp$r - 30))]
  }
  else{
    lung_meta$k_Other[i] <- NA
  }

  # # CD 19
  if(lung_meta$p_CD19[i] > 0){
    # k function estimated for image i B cells using isotropic correction
    kEstCD19_temp <- MItools::permute.envelope(pppObj, fun = Kest, funargs = list(i = "CD19", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of 30 microns
    lung_meta$k_CD19[i] <- (kEstCD19_temp$obs - kEstCD19_temp$mmean)[
      which.min(abs(kEstCD19_temp$r - 30))]
  }
  else{
    lung_meta$k_CD19[i] <- NA
  }

  # # CD4
  if(lung_meta$p_CD4[i] > 0){
    kEstCD4_temp <- permute.envelope(pppObj, fun = Kest, funargs = list(i = "CD4", correction = "isotropic"))
    lung_meta$k_CD4[i] <- (kEstCD4_temp$obs - kEstCD4_temp$mmean)[
      which.min(abs(kEstCD4_temp$r - 30))]
  }
  else{
    lung_meta$k_CD14[i] <- NA
  }
  
  # # Calculating K-cross values
  # # CK x CD8
  if(lung_meta$p_CD8[i] > 0 & lung_meta$p_CK[i] > 0){
    kCross_temp <- permute.envelope(pppObj, fun = Kcross, funargs = list(i = "CK", j = "CD8", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of 30 microns
    lung_meta$k_CK_CD8[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 30))]
  }
  else{
    lung_meta$k_CK_CD8[i] <- NA
  }

  # CK x CD14
  if(lung_meta$p_CD14[i] > 0 & lung_meta$p_CK[i] > 0){
    kCross_temp <- permute.envelope(pppObj, fun = Kcross, funargs = list(i = "CK", j = "CD14", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of 30 microns
    lung_meta$k_CK_CD14[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 30))]
  }
  else{
    lung_meta$k_CK_CD14[i] <- NA
  }

  # CK x Other
  if(lung_meta$p_Other[i] > 0 & lung_meta$p_CK[i] > 0){
    kCross_temp <- permute.envelope(pppObj, fun = Kcross, funargs = list(i = "CK", j = "Other", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of 30 microns
    lung_meta$k_CK_Other[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 30))]
  }
  else{
    lung_meta$k_CK_Other[i] <- NA
  }
  
  # CK x CD19
  if(lung_meta$p_CD19[i] > 0 & lung_meta$p_CK[i] > 0){
    # k function estimated for image i B cells using isotropic correction
    kCross_temp <- MItools::permute.envelope(pppObj, fun = Kcross, funargs = list(i = "CK", j = "CD19", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of 30 microns
    lung_meta$k_CK_CD19[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 30))
    ]
  }
  else{
    # if there are no CD19 cells or no tumor cells, record cross function value as NA
    lung_meta$k_CK_CD19[i] <- NA
  }

  # CK x CD4
  if(lung_meta$p_CD4[i] > 0 & lung_meta$p_CK[i] > 0){
    kCross_temp <- permute.envelope(pppObj, fun = Kcross, funargs = list(i = "CK", j = "CD4", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of 30 microns
    lung_meta$k_CK_CD4[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 30))]
  }
  else{
    lung_meta$k_CK_CD4[i] <- NA
  }

  # CD8 x CD14
  if(lung_meta$p_CD8[i] > 0 & lung_meta$p_CD14[i] > 0){
    kCross_temp <- permute.envelope(pppObj, fun = Kcross, funargs = list(i = "CD14", j = "CD8", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of 30 microns
    lung_meta$k_CD8_CD14[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 30))]
  }
  else{
    lung_meta$k_CD8_CD14[i] <- NA
  }

  # CD8 x Other
  if(lung_meta$p_CD8[i] > 0 & lung_meta$p_Other[i] > 0){
    kCross_temp <- permute.envelope(pppObj, fun = Kcross, funargs = list(i = "Other", j = "CD8", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of 30 microns
    lung_meta$k_CD8_Other[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 30))]
  }
  else{
    lung_meta$k_CD8_Other[i] <- NA
  }

  # CD8 x CD19
  if(lung_meta$p_CD8[i] > 0 & lung_meta$p_CD19[i] > 0){
    kCross_temp <- permute.envelope(pppObj, fun = Kcross, funargs = list(i = "CD19", j = "CD8", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of 30 microns
    lung_meta$k_CD8_CD19[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 30))]
  }
  else{
    lung_meta$k_CD8_CD19[i] <- NA
  }

  # CD8 x CD4
  if(lung_meta$p_CD8[i] > 0 & lung_meta$p_CD4[i] > 0){
    kCross_temp <- permute.envelope(pppObj, fun = Kcross, funargs = list(i = "CD4", j = "CD8", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of 30 microns
    lung_meta$k_CD8_CD4[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 30))]
  }
  else{
    lung_meta$k_CD8_CD4[i] <- NA
  }

  # CD14 x Other
  if(lung_meta$p_Other[i] > 0 & lung_meta$p_CD14[i] > 0){
    kCross_temp <- permute.envelope(pppObj, fun = Kcross, funargs = list(i = "CD14", j = "Other", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of 30 microns
    lung_meta$k_CD14_Other[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 30))]
  }
  else{
    lung_meta$k_CD14_Other[i] <- NA
  }

  # CD14 x CD19
  if(lung_meta$p_CD19[i] > 0 & lung_meta$p_CD14[i] > 0){
    kCross_temp <- permute.envelope(pppObj, fun = Kcross, funargs = list(i = "CD14", j = "CD19", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of 30 microns
    lung_meta$k_CD14_CD19[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 30))]
  }
  else{
    lung_meta$k_CD14_CD19[i] <- NA
  }

  # CD14 x CD4
  if(lung_meta$p_CD4[i] > 0 & lung_meta$p_CD14[i] > 0){
    kCross_temp <- permute.envelope(pppObj, fun = Kcross, funargs = list(i = "CD14", j = "CD4", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of 30 microns
    lung_meta$k_CD14_CD4[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 30))]
  }
  else{
    lung_meta$k_CD14_CD4[i] <- NA
  }

  # Other x CD19
  if(lung_meta$p_CD19[i] > 0 & lung_meta$p_Other[i] > 0){
    kCross_temp <- permute.envelope(pppObj, fun = Kcross, funargs = list(i = "CD19", j = "Other", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of 30 microns
    lung_meta$k_Other_CD19[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 30))]
  }
  else{
    lung_meta$k_Other_CD19[i] <- NA
  }

  # Other x CD4
  if(lung_meta$p_CD4[i] > 0 & lung_meta$p_Other[i] > 0){
    kCross_temp <- permute.envelope(pppObj, fun = Kcross, funargs = list(i = "CD4", j = "Other", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of 30 microns
    lung_meta$k_Other_CD4[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 30))]
  }
  else{
    lung_meta$k_Other_CD4[i] <- NA
  }

  # CD19 x CD4
  if(lung_meta$p_CD4[i] > 0 & lung_meta$p_CD19[i] > 0){
    kCross_temp <- permute.envelope(pppObj, fun = Kcross, funargs = list(i = "CD4", j = "CD19", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of 30 microns
    lung_meta$k_CD19_CD4[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - 30))]
  }
  else{
    lung_meta$k_CD19_CD4[i] <- NA
  }
}

# filtering to only include images with most number of cells per patient
lung_meta30 <- lung_meta |>
  filter(image_id %in% lung_filtered$image_id)
  
