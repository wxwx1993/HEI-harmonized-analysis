library("gnm")

dir_data = "/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/"

load(paste0(dir_data,"aggregate_data_rm.RData"))

# cut-off
# < 5 
aggregate_data_rm_sub5 <- subset(aggregate_data_rm, pm25 < 5)
gnm_pm25_mod4.5 <- gnm(dead ~ pm25 + as.factor(year) +
                       medhouseholdincome + medianhousevalue + poverty + education + pct_owner_occ +
                       as.factor(region) + 
                       offset(log(time_count)),
                     eliminate = (as.factor(sex):as.factor(entry_age_break):as.factor(race):as.factor(dual):as.factor(followup_year)),
                     data = aggregate_data_rm_sub5,
                     family = poisson(link = "log"))
summary_pm25_mod4.5 <- summary(gnm_pm25_mod4.5)
rm(gnm_pm25_mod4.5)
exp(summary_pm25_mod4.5$coefficients[1]*10)
#[1] 1.051099
save(summary_pm25_mod4.5, file="/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/summary_pm25_mod4.5.RData")

# < 7.5 
aggregate_data_rm_sub7.5 <- subset(aggregate_data_rm, pm25 < 7.5)
gnm_pm25_mod4.7.5 <- gnm(dead ~ pm25 + as.factor(year) +
                         medhouseholdincome + medianhousevalue + poverty + education + pct_owner_occ +
                         as.factor(region) + 
                         offset(log(time_count)),
                       eliminate = (as.factor(sex):as.factor(entry_age_break):as.factor(race):as.factor(dual):as.factor(followup_year)),
                       data = aggregate_data_rm_sub7.5,
                       family = poisson(link = "log"))
summary_pm25_mod4.7.5 <- summary(gnm_pm25_mod4.7.5)
rm(gnm_pm25_mod4.7.5)
exp(summary_pm25_mod4.7.5$coefficients[1]*10)
#[1] 
save(summary_pm25_mod4.7.5, file="/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/summary_pm25_mod4.7.5.RData")

# < 10 
aggregate_data_rm_sub10 <- subset(aggregate_data_rm, pm25 < 10)
gnm_pm25_mod4.10 <- gnm(dead ~ pm25 + as.factor(year) +
                           medhouseholdincome + medianhousevalue + poverty + education + pct_owner_occ +
                           as.factor(region) + 
                           offset(log(time_count)),
                         eliminate = (as.factor(sex):as.factor(entry_age_break):as.factor(race):as.factor(dual):as.factor(followup_year)),
                         data = aggregate_data_rm_sub10,
                         family = poisson(link = "log"))
summary_pm25_mod4.10 <- summary(gnm_pm25_mod4.10)
rm(gnm_pm25_mod4.10)
exp(summary_pm25_mod4.10$coefficients[1]*10)
#[1] 
save(summary_pm25_mod4.10, file="/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/summary_pm25_mod4.10.RData")


# < 12 
aggregate_data_rm_sub12 <- subset(aggregate_data_rm, pm25 < 12)
gnm_pm25_mod4.12 <- gnm(dead ~ pm25 + as.factor(year) +
                           medhouseholdincome + medianhousevalue + poverty + education + pct_owner_occ +
                           as.factor(region) + 
                           offset(log(time_count)),
                         eliminate = (as.factor(sex):as.factor(entry_age_break):as.factor(race):as.factor(dual):as.factor(followup_year)),
                         data = aggregate_data_rm12,
                         family = poisson(link = "log"))
summary_pm25_mod4.12 <- summary(gnm_pm25_mod4.12)
rm(gnm_pm25_mod4.12)
exp(summary_pm25_mod4.12$coefficients[1]*10)
#[1] 
save(summary_pm25_mod4.12, file="/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/summary_pm25_mod4.12.RData")

