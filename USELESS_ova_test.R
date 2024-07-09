
new_ovarian_meta40 = subset(ovarian_meta40, select = -c(sample_id, tma, diagnosis, primary, recurrent, treatment_effect, stage, grade, survival_time, death, BRCA_mutation, age_at_diagnosis, time_to_recurrence, parpi_inhibitor, debulking, totalCell, p_Other, p_Tumor, p_BCell, p_Macrophage, p_THelper, p_CytotoxicT))

new_ovarian_meta40 = na.omit(new_ovarian_meta40)

new_ovarian_meta40 = scale(new_ovarian_meta40)

View(new_ovarian_meta40)


View(ovarian_meta40)