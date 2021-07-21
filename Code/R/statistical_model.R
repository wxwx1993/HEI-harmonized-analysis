library("gnm")

dir_data = "/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/"

load(paste0(dir_data,"aggregate_data_rm.RData"))

#### model1
gnm_pm25_mod1 <- gnm(dead ~ pm25 + as.factor(year) +
                       offset(log(time_count)),
                     eliminate = (as.factor(sex):as.factor(entry_age_break):as.factor(followup_year)),
                     data = aggregate_data_rm,
                     family = poisson(link = "log"))
summary_pm25_mod1 <- summary(gnm_pm25_mod1)
rm(gnm_pm25_mod1)
exp(summary_pm25_mod1$coefficients[1]*10)
#[1] 1.035655

save(summary_pm25_mod1, file="/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/summary_pm25_mod1.RData")

# alternative specifications (unused)
gnm_pm25_mod1.1 <- gnm(dead ~ pm25 + as.factor(year) + as.factor(sex) + 
                       offset(log(time_count)),
                     eliminate = (as.factor(entry_age_break):as.factor(followup_year)),
                     data = aggregate_data_rm,
                     family = poisson(link = "log"))
summary_pm25_mod1.1 <- summary(gnm_pm25_mod1.1)
rm(gnm_pm25_mod1.1)
save(summary_pm25_mod1.1, file="/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/summary_pm25_mod1.1.RData")

gnm_pm25_mod1.2 <- gnm(dead ~ pm25 + as.factor(year) + as.factor(sex) + as.factor(entry_age_break) + 
                         offset(log(time_count)),
                       eliminate = as.factor(followup_year),
                       data = aggregate_data_rm,
                       family = poisson(link = "log"))
summary_pm25_mod1.2 <- summary(gnm_pm25_mod1.2)
rm(gnm_pm25_mod1.2)
save(summary_pm25_mod1.2, file="/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/summary_pm25_mod1.2.RData")

#### model2
gnm_pm25_mod2 <- gnm(dead ~ pm25 + as.factor(year) +
                       offset(log(time_count)),
                     eliminate = (as.factor(sex):as.factor(entry_age_break):as.factor(race):as.factor(dual):as.factor(followup_year)),
                     data = aggregate_data_rm,
                     family = poisson(link = "log"))
summary_pm25_mod2 <- summary(gnm_pm25_mod2)
rm(gnm_pm25_mod2)
exp(summary_pm25_mod2$coefficients[1]*10)
#[1] 1.025955
save(summary_pm25_mod2, file="/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/summary_pm25_mod2.RData")

#### model3
gnm_pm25_mod3 <- gnm(dead ~ pm25 + as.factor(year) + 
                       medhouseholdincome + medianhousevalue + poverty + education + pct_owner_occ +
                       offset(log(time_count)),
                     eliminate = (as.factor(sex):as.factor(entry_age_break):as.factor(race):as.factor(dual):as.factor(followup_year)),
                     data = aggregate_data_rm,
                     family = poisson(link = "log"))
summary_pm25_mod3 <- summary(gnm_pm25_mod3)
rm(gnm_pm25_mod3)
exp(summary_pm25_mod3$coefficients[1]*10)
#[1] 1.054755
save(summary_pm25_mod3, file="/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/summary_pm25_mod3.RData")

#### model4
gnm_pm25_mod4 <- gnm(dead ~ pm25 + as.factor(year) +
                       medhouseholdincome + medianhousevalue + poverty + education + pct_owner_occ +
                       as.factor(region) + 
                       offset(log(time_count)),
                     eliminate = (as.factor(sex):as.factor(entry_age_break):as.factor(race):as.factor(dual):as.factor(followup_year)),
                     data = aggregate_data_rm,
                     family = poisson(link = "log"))
summary_pm25_mod4 <- summary(gnm_pm25_mod4)
rm(gnm_pm25_mod4)
exp(summary_pm25_mod4$coefficients[1]*10)
#[1] 1.051028
save(summary_pm25_mod4, file="/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/summary_pm25_mod4.RData")


