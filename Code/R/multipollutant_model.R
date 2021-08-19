library("gnm")
library("stringr")
library("mgcv")
library("dplyr")

dir_data = "/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/"

load(paste0(dir_data, "aggregate_data_rm.RData"))

allyears_ozone <- read.csv("/nfs/home/X/xwu/shared_space/ci3_exposure/ozone/whole_us/annual/zipcode/requaia_predictions/ywei_aggregation/all_years.csv")
allyears_ozone$ZIP <- str_pad(allyears_ozone$ZIP, 5, pad = "0")
allyears_no2 <- read.csv("/nfs/home/X/xwu/shared_space/ci3_exposure/no2/whole_us/annual/zipcode/qd_predictions_ensemble/ywei_aggregations/all_years.csv")
allyears_no2$ZIP <- str_pad(allyears_no2$ZIP, 5, pad = "0")

aggregate_data_rm = merge(aggregate_data_rm, 
                          allyears_ozone, 
                          by.x = c("zip","year"), 
                          by.y = c("ZIP","year"))
aggregate_data_rm = merge(aggregate_data_rm, 
                          allyears_no2, 
                          by.x = c("zip","year"), 
                          by.y = c("ZIP","year"))

aggregate_data_rm$ox = (1.07 * aggregate_data_rm$no2 + 2.075 * aggregate_data_rm$ozone) / 3.145

## Classic two-pollutant
# pm + no2
gnm_pm25_mod4_no2 <- gnm(dead ~ pm25 + no2 + as.factor(year) +
                         medhouseholdincome + medianhousevalue + poverty + education + pct_owner_occ +
                         as.factor(region) + 
                         offset(log(time_count)),
                       eliminate = (as.factor(sex):as.factor(entry_age_break):as.factor(race):as.factor(dual):as.factor(followup_year)),
                       data = aggregate_data_rm,
                       family = poisson(link = "log"))
summary_pm25_mod4_no2 <- summary(gnm_pm25_mod4_no2)
rm(gnm_pm25_mod4_no2)
exp(summary_pm25_mod4_no2$coefficients[1:2]*10)
#[1] 
save(summary_pm25_mod4_no2, file="/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/summary_pm25_mod4_no2.RData")

# pm + o3
gnm_pm25_mod4_ozone <- gnm(dead ~ pm25 + ozone + as.factor(year) +
                           medhouseholdincome + medianhousevalue + poverty + education + pct_owner_occ +
                           as.factor(region) + 
                           offset(log(time_count)),
                         eliminate = (as.factor(sex):as.factor(entry_age_break):as.factor(race):as.factor(dual):as.factor(followup_year)),
                         data = aggregate_data_rm,
                         family = poisson(link = "log"))
summary_pm25_mod4_ozone <- summary(gnm_pm25_mod4_ozone)
rm(gnm_pm25_mod4_ozone)
exp(summary_pm25_mod4.7.5$coefficients[1:2]*10)
#[1] 
save(summary_pm25_mod4_ozone, file="/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/summary_pm25_mod4_ozone.RData")

# pm + ox
gnm_pm25_mod4_ox <- gnm(dead ~ pm25 + ox + as.factor(year) +
                          medhouseholdincome + medianhousevalue + poverty + education + pct_owner_occ +
                          as.factor(region) + 
                          offset(log(time_count)),
                        eliminate = (as.factor(sex):as.factor(entry_age_break):as.factor(race):as.factor(dual):as.factor(followup_year)),
                        data = aggregate_data_rm,
                        family = poisson(link = "log"))
summary_pm25_mod4_ox <- summary(gnm_pm25_mod4_ox)
rm(gnm_pm25_mod4_ox)
exp(summary_pm25_mod4.10$coefficients[1:2]*10)
#[1] 
save(summary_pm25_mod4_ox, file="/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/summary_pm25_mod4_ox.RData")


## Bivariate thin plate spline surface
# pm25 + no2
gam_pm25_mod4_no2 <- gam(dead ~ s(pm25, no2) + as.factor(year) +
                          medhouseholdincome + medianhousevalue + poverty + education + pct_owner_occ +
                          as.factor(region) + 
                          offset(log(time_count)) +
                         (as.factor(sex):as.factor(entry_age_break):as.factor(race):as.factor(dual):as.factor(followup_year)),
                        data = aggregate_data_rm,
                        family = poisson(link = "log"))
plot(gam_pm25_mod4_no2)
#[1] 
save(gam_pm25_mod4_no2, file="/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/gam_pm25_mod4_no2.RData")

# pm25 + ozone
gam_pm25_mod4_ozone <- gam(dead ~ s(pm25, ozone) + as.factor(year) +
                           medhouseholdincome + medianhousevalue + poverty + education + pct_owner_occ +
                           as.factor(region) + 
                           offset(log(time_count)) +
                           (as.factor(sex):as.factor(entry_age_break):as.factor(race):as.factor(dual):as.factor(followup_year)),
                         data = aggregate_data_rm,
                         family = poisson(link = "log"))
plot(gam_pm25_mod4_ozone)
#[1] 
save(gam_pm25_mod4_ozone, file="/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/gam_pm25_mod4_ozone.RData")

# pm25 + ox
gam_pm25_mod4_ox <- gam(dead ~ s(pm25, ox) + as.factor(year) +
                             medhouseholdincome + medianhousevalue + poverty + education + pct_owner_occ +
                             as.factor(region) + 
                             offset(log(time_count)) +
                             (as.factor(sex):as.factor(entry_age_break):as.factor(race):as.factor(dual):as.factor(followup_year)),
                           data = aggregate_data_rm,
                           family = poisson(link = "log"))
plot(gam_pm25_mod4_ox)
#[1] 
save(gam_pm25_mod4_ox, file="/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/gam_pm25_mod4_ox.RData")

# Examine association between PM2.5 and mortality by tertile of the other pollutant
# align with MAPLE team to use tertiles by person-year exposure
# tertile by ozone (the following code could be optimized)
person_year_ozone <- aggregate_data_rm %>% 
  group_by(zip, year) %>% slice(c(1)) 
q <- quantile(person_year_ozone$ozone, c(1/3,2/3))
person_year_ozone$tertiles <- NA
person_year_ozone$tertiles[person_year_ozone$ozone <= q[1]] <- 1
person_year_ozone$tertiles[person_year_ozone$ozone > q[1] & person_year_ozone$ozone <= q[2]] <- 2
person_year_ozone$tertiles[person_year_ozone$ozone > q[2]] <- 3

aggregate_data_rm <- merge(aggregate_data_rm, 
                           person_year_ozone[,c("zip", "year", "tertiles")], 
                           by.x = c("zip","year"), 
                           by.y = c("zip","year"))

gnm_pm25_mod4_ozone_tert1 <- gnm(dead ~ pm25 + as.factor(year) +
                           medhouseholdincome + medianhousevalue + poverty + education + pct_owner_occ +
                           as.factor(region) + 
                           offset(log(time_count)),
                         eliminate = (as.factor(sex):as.factor(entry_age_break):as.factor(race):as.factor(dual):as.factor(followup_year)),
                         data = subset(aggregate_data_rm,tertiles == 1),
                         family = poisson(link = "log"))
summary_gnm_pm25_mod4_ozone_tert1 <- summary(gnm_pm25_mod4_ozone_tert1)
rm(gnm_pm25_mod4_ozone_tert1)
exp(summary_gnm_pm25_mod4_ozone_tert1$coefficients[1]*10)
#
save(summary_gnm_pm25_mod4_ozone_tert1, file="/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/summary_gnm_pm25_mod4_ozone_tert1.RData")


gnm_pm25_mod4_ozone_tert2 <- gnm(dead ~ pm25 + as.factor(year) +
                                   medhouseholdincome + medianhousevalue + poverty + education + pct_owner_occ +
                                   as.factor(region) + 
                                   offset(log(time_count)),
                                 eliminate = (as.factor(sex):as.factor(entry_age_break):as.factor(race):as.factor(dual):as.factor(followup_year)),
                                 data = subset(aggregate_data_rm,tertiles == 2),
                                 family = poisson(link = "log"))
summary_gnm_pm25_mod4_ozone_tert2 <- summary(gnm_pm25_mod4_ozone_tert2)
rm(gnm_pm25_mod4_ozone_tert2)
exp(summary_gnm_pm25_mod4_ozone_tert2$coefficients[1]*10)
#
save(summary_gnm_pm25_mod4_ozone_tert2, file="/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/summary_gnm_pm25_mod4_ozone_tert2.RData")

gnm_pm25_mod4_ozone_tert3 <- gnm(dead ~ pm25 + as.factor(year) +
                                   medhouseholdincome + medianhousevalue + poverty + education + pct_owner_occ +
                                   as.factor(region) + 
                                   offset(log(time_count)),
                                 eliminate = (as.factor(sex):as.factor(entry_age_break):as.factor(race):as.factor(dual):as.factor(followup_year)),
                                 data = subset(aggregate_data_rm,tertiles == 3),
                                 family = poisson(link = "log"))
summary_gnm_pm25_mod4_ozone_tert3 <- summary(gnm_pm25_mod4_ozone_tert3)
rm(gnm_pm25_mod4_ozone_tert3)
exp(summary_gnm_pm25_mod4_ozone_tert3$coefficients[1]*10)
#
save(summary_gnm_pm25_mod4_ozone_tert3, file="/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/summary_gnm_pm25_mod4_ozone_tert3.RData")

# tertile by no2
# TBD (use the same code as the previous chunk but switch ozone to no2)


# tertile by ox
# TBD (use the same code as the previous chunk but switch ozone to ox)


