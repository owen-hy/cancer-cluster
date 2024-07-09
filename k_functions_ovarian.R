library(tidyverse)
#devtools::install_github("junsoukchoi/MItools")
library(MItools)

ovarianData <- read_rds("./ovarianData.RDS")
ovarian_cell <- ovarianData$ovarian_cells
ovarian_meta <- ovarianData$ovarian_metadata

# The micron value of R we use to calculate our K/K Cross function values 
micron <- 30

cell_types <- unique(ovarian_cell$pheno) # 6 unqiue cell types (Other, Tumor, B Cell, Macrophage, 
                                         # T Helper, Cytotoxic T)

# add image-level prevalence of cell types to metadata
ovarian_meta <- ovarian_meta |>
  left_join(ovarian_cell |>
              group_by(sample_id) |>
              summarize(totalCell = n(),
                        p_Other = sum(pheno == "Other")/n(),
                        p_Tumor = sum(pheno == "Tumor")/n(),
                        p_BCell = sum(pheno == "B Cell")/n(),
                        p_Macrophage = sum(pheno == "Macrophage")/n(),
                        p_THelper = sum(pheno == "T Helper")/n(),
                        p_CytotoxicT = sum(pheno == "Cytotoxic T")/n()),
            join_by(sample_id))


# create variables for k function values of each image
ovarian_meta$k_Other <- NA # k function values for Other
ovarian_meta$k_Tumor <- NA # k function values for Tumor
ovarian_meta$k_BCell <- NA # k function values for B Cell
ovarian_meta$k_Macrophage <- NA # k function values for Macrophage
ovarian_meta$k_THelper <- NA # k function values for T Helper
ovarian_meta$k_CytotoxicT <- NA # k function values for Cytotoxic T

# create variables for k cross values of each image
ovarian_meta$k_Other_Tumor <- NA # k cross values for Other and Tumor
ovarian_meta$k_Other_BCell <- NA # k cross values for Other and B Cell
ovarian_meta$k_Other_Macrophage <- NA # k cross values for Other and Macrophage
ovarian_meta$k_Other_THelper <- NA # k cross values for Other and T Helper
ovarian_meta$k_Other_CytotoxicT <- NA # k cross values for Other and Cytotoxic T
ovarian_meta$k_Tumor_BCell <- NA # k cross values for Tumor and B Cell
ovarian_meta$k_Tumor_Macrophage <- NA # k cross values for Tumor and Macrophage
ovarian_meta$k_Tumor_THelper <- NA # k cross values for Tumor and T Helper
ovarian_meta$k_Tumor_CytotoxicT <- NA # k cross values for Tumor and Cytotoxic T
ovarian_meta$k_BCell_Macrophage <- NA # k cross values for B Cell and Macrophage
ovarian_meta$k_BCell_THelper <- NA # k cross values for B Cell and T Helper
ovarian_meta$k_BCell_CytotoxicT <- NA # k cross values for B Cell and Cytotoxic T
ovarian_meta$k_Macrophage_THelper <- NA # k cross values for Macrophage and T Helper
ovarian_meta$k_Macrophage_CytotoxicT <- NA # k cross values for Macrophage and Cytotoxic T
ovarian_meta$k_THelper_CytotoxicT <- NA # k cross values for T Helper and Cytotoxic T

# for loop to iterate through each image in the metadata
for(i in 1:nrow(ovarian_meta)){
  # temporary dataset for image "i" in the ith iteration of this "for" loop
  tempData <- ovarianData$ovarian_cells |>
    filter(sample_id == ovarian_meta$sample_id[i])
  
  # point pattern object for all of image i
  pppObj <- spatstat.geom::ppp(
    x = tempData$x, y = tempData$y,
    xrange = range(tempData$x), yrange = range(tempData$y),
    marks = factor(tempData$pheno))
  
  if(ovarian_meta$p_Other[i] > 0){
    # k function estimated for image i Misc. cells using isotropic correction
    kEstOther_temp <- permute.envelope(pppObj, fun = Kest, funargs = list(i = "Other", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of micron microns
    ovarian_meta$k_Other[i] <- (kEstOther_temp$obs - kEstOther_temp$mmean)[
      which.min(abs(kEstOther_temp$r - micron))]
  }
  else{
    # if there are no tumor cells, record NA for this image's k function
    ovarian_meta$k_Other[i] <- NA
  }
  
  # Tumor
  if(ovarian_meta$p_Tumor[i] > 0){
    kEstTumor_temp <- permute.envelope(pppObj, fun = Kest, funargs = list(i = "Tumor", correction = "isotropic"))
    ovarian_meta$k_Tumor[i] <- (kEstTumor_temp$obs - kEstTumor_temp$mmean)[
      which.min(abs(kEstTumor_temp$r - micron))]
  }
  else{
    ovarian_meta$k_Tumor[i] <- NA
  }
  
  # # B Cell
  if(ovarian_meta$p_BCell[i] > 0){
    kEstBCell_temp <- permute.envelope(pppObj, fun = Kest, funargs = list(i = "B Cell", correction = "isotropic"))
    ovarian_meta$k_BCell[i] <- (kEstBCell_temp$obs - kEstBCell_temp$mmean)[
      which.min(abs(kEstBCell_temp$r - micron))]
  }
  else{
    ovarian_meta$k_BCell[i] <- NA
  }
  
  # # Macrophage
  if(ovarian_meta$p_Macrophage[i] > 0){
    kEstMacrophage_temp <- permute.envelope(pppObj, fun = Kest, funargs = list(i = "Macrophage", correction = "isotropic"))
    ovarian_meta$k_Macrophage[i] <- (kEstMacrophage_temp$obs - kEstMacrophage_temp$mmean)[
      which.min(abs(kEstMacrophage_temp$r - micron))]
  }
  else{
    ovarian_meta$k_Macrophage[i] <- NA
  }
  
  # # T Helper
  if(ovarian_meta$p_THelper[i] > 0){
    # k function estimated for image i B cells using isotropic correction
    kEstTHelper_temp <- MItools::permute.envelope(pppObj, fun = Kest, funargs = list(i = "T Helper", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of micron microns
    ovarian_meta$k_THelper[i] <- (kEstTHelper_temp$obs - kEstTHelper_temp$mmean)[
      which.min(abs(kEstTHelper_temp$r - micron))]
  }
  else{
    ovarian_meta$k_THelper[i] <- NA
  }
  
  # # Cytotoxic T
  if(ovarian_meta$p_CytotoxicT[i] > 0){
    kEstCytotoxicT_temp <- permute.envelope(pppObj, fun = Kest, funargs = list(i = "Cytotoxic T", correction = "isotropic"))
    ovarian_meta$k_CytotoxicT[i] <- (kEstCytotoxicT_temp$obs - kEstCytotoxicT_temp$mmean)[
      which.min(abs(kEstCytotoxicT_temp$r - micron))]
  }
  else{
    ovarian_meta$k_CytotoxicT[i] <- NA
  }
  
  # # Calculating K-cross values
  # # Other x Tumor
  if(ovarian_meta$p_Tumor[i] > 0 & ovarian_meta$p_Other[i] > 0){
    kCross_temp <- permute.envelope(pppObj, fun = Kcross, funargs = list(i = "Other", j = "Tumor", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of micron microns
    ovarian_meta$k_Other_Tumor[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - micron))]
  }
  else{
    ovarian_meta$k_Other_Tumor[i] <- NA
  }
  
  # Other x B Cell
  if(ovarian_meta$p_BCell[i] > 0 & ovarian_meta$p_Other[i] > 0){
    kCross_temp <- permute.envelope(pppObj, fun = Kcross, funargs = list(i = "Other", j = "B Cell", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of micron microns
    ovarian_meta$k_Other_BCell[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - micron))]
  }
  else{
    ovarian_meta$k_Other_BCell[i] <- NA
  }
  
  # Other x Macrophage
  if(ovarian_meta$p_Other[i] > 0 & ovarian_meta$p_Macrophage[i] > 0){
    kCross_temp <- permute.envelope(pppObj, fun = Kcross, funargs = list(i = "Other", j = "Macrophage", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of micron microns
    ovarian_meta$k_Other_Macrophage[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - micron))]
  }
  else{
    ovarian_meta$k_Other_Macrophage[i] <- NA
  }
  
  # Other x THelper
  if(ovarian_meta$p_THelper[i] > 0 & ovarian_meta$p_Other[i] > 0){
    # k function estimated for image i B cells using isotropic correction
    kCross_temp <- MItools::permute.envelope(pppObj, fun = Kcross, funargs = list(i = "Other", j = "T Helper", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of micron microns
    ovarian_meta$k_Other_THelper[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - micron))
    ]
  }
  else{
    # if there are no T Helper cells or no tumor cells, record cross function value as NA
    ovarian_meta$k_Other_THelper[i] <- NA
  }
  
  # Other x CytotoxicT
  if(ovarian_meta$p_CytotoxicT[i] > 0 & ovarian_meta$p_Other[i] > 0){
    kCross_temp <- permute.envelope(pppObj, fun = Kcross, funargs = list(i = "Other", j = "Cytotoxic T", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of micron microns
    ovarian_meta$k_Other_CytotoxicT[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - micron))]
  }
  else{
    ovarian_meta$k_Other_CytotoxicT[i] <- NA
  }
  
  # Tumor x B Cell
  if(ovarian_meta$p_Tumor[i] > 0 & ovarian_meta$p_BCell[i] > 0){
    kCross_temp <- permute.envelope(pppObj, fun = Kcross, funargs = list(i = "B Cell", j = "Tumor", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of micron microns
    ovarian_meta$k_Tumor_BCell[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - micron))]
  }
  else{
    ovarian_meta$k_Tumor_BCell[i] <- NA
  }
  
  # Tumor x Macrophage
  if(ovarian_meta$p_Tumor[i] > 0 & ovarian_meta$p_Macrophage[i] > 0){
    kCross_temp <- permute.envelope(pppObj, fun = Kcross, funargs = list(i = "Macrophage", j = "Tumor", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of micron microns
    ovarian_meta$k_Tumor_Macrophage[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - micron))]
  }
  else{
    ovarian_meta$k_Tumor_Macrophage[i] <- NA
  }
  
  # Tumor x THelper
  if(ovarian_meta$p_Tumor[i] > 0 & ovarian_meta$p_THelper[i] > 0){
    kCross_temp <- permute.envelope(pppObj, fun = Kcross, funargs = list(i = "T Helper", j = "Tumor", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of micron microns
    ovarian_meta$k_Tumor_THelper[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - micron))]
  }
  else{
    ovarian_meta$k_Tumor_THelper[i] <- NA
  }
  
  # Tumor x CytotoxicT
  if(ovarian_meta$p_Tumor[i] > 0 & ovarian_meta$p_CytotoxicT[i] > 0){
    kCross_temp <- permute.envelope(pppObj, fun = Kcross, funargs = list(i = "Cytotoxic T", j = "Tumor", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of micron microns
    ovarian_meta$k_Tumor_CytotoxicT[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - micron))]
  }
  else{
    ovarian_meta$k_Tumor_CytotoxicT[i] <- NA
  }
  
  # B Cell x Macrophage
  if(ovarian_meta$p_Macrophage[i] > 0 & ovarian_meta$p_BCell[i] > 0){
    kCross_temp <- permute.envelope(pppObj, fun = Kcross, funargs = list(i = "B Cell", j = "Macrophage", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of micron microns
    ovarian_meta$k_BCell_Macrophage[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - micron))]
  }
  else{
    ovarian_meta$k_BCell_Macrophage[i] <- NA
  }
  
  # B Cell x THelper
  if(ovarian_meta$p_THelper[i] > 0 & ovarian_meta$p_BCell[i] > 0){
    kCross_temp <- permute.envelope(pppObj, fun = Kcross, funargs = list(i = "B Cell", j = "T Helper", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of micron microns
    ovarian_meta$k_BCell_THelper[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - micron))]
  }
  else{
    ovarian_meta$k_BCell_THelper[i] <- NA
  }
  
  # B Cell x CytotoxicT
  if(ovarian_meta$p_CytotoxicT[i] > 0 & ovarian_meta$p_BCell[i] > 0){
    kCross_temp <- permute.envelope(pppObj, fun = Kcross, funargs = list(i = "B Cell", j = "Cytotoxic T", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of micron microns
    ovarian_meta$k_BCell_CytotoxicT[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - micron))]
  }
  else{
    ovarian_meta$k_BCell_CytotoxicT[i] <- NA
  }
  
  # Macrophage x THelper
  if(ovarian_meta$p_THelper[i] > 0 & ovarian_meta$p_Macrophage[i] > 0){
    kCross_temp <- permute.envelope(pppObj, fun = Kcross, funargs = list(i = "T Helper", j = "Macrophage", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of micron microns
    ovarian_meta$k_Macrophage_THelper[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - micron))]
  }
  else{
    ovarian_meta$k_Macrophage_THelper[i] <- NA
  }
  
  # Macrophage x CytotoxicT
  if(ovarian_meta$p_CytotoxicT[i] > 0 & ovarian_meta$p_Macrophage[i] > 0){
    kCross_temp <- permute.envelope(pppObj, fun = Kcross, funargs = list(i = "Cytotoxic T", j = "Macrophage", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of micron microns
    ovarian_meta$k_Macrophage_CytotoxicT[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - micron))]
  }
  else{
    ovarian_meta$k_Macrophage_CytotoxicT[i] <- NA
  }
  
  # THelper x CytotoxicT
  if(ovarian_meta$p_CytotoxicT[i] > 0 & ovarian_meta$p_THelper[i] > 0){
    kCross_temp <- permute.envelope(pppObj, fun = Kcross, funargs = list(i = "Cytotoxic T", j = "T Helper", correction = "isotropic"))
    # get difference between estimate and theoretical value at radius of micron microns
    ovarian_meta$k_THelper_CytotoxicT[i] <- (kCross_temp$obs - kCross_temp$mmean)[
      which.min(abs(kCross_temp$r - micron))]
  }
  else{
    ovarian_meta$k_THelper_CytotoxicT[i] <- NA
  }
}


