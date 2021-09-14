library("gnm")
library("stringr")
library("mgcv")
library("dplyr")
library("CausalGPS")

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

aggregate_data_rm$year<-as.factor(aggregate_data_rm$year)
aggregate_data_rm$region<-as.factor(aggregate_data_rm$region)
aggregate_data_rm$zip<-as.character(aggregate_data_rm$zip)


zip_year <- aggregate_data_rm %>% 
  group_by(zip, year) %>% dplyr::slice(c(1)) 

a.vals <- seq(min(zip_year$pm25), max(zip_year$pm25), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])

# trim 0.1/0.99 of the matched set for stabilty (error if no trim)
# GPS: PM | potential confounders
match_pop_all <- generate_pseudo_pop(Y = zip_year$zip,
                                     w = zip_year$pm25,
                                     c = zip_year[, c(2, 10:13,15:16)],
                                     ci_appr = "matching",
                                     pred_model = "sl",
                                     gps_model = "parametric",
                                     use_cov_transform = FALSE,
                                     #transformers = list("pow2", "pow3"),
                                     sl_lib = c("m_xgboost"),
                                     params = list("xgb_nrounds" = 50,
                                                   "xgb_max_depth" = 6,
                                                   "xgb_eta" = 0.3,
                                                   "xgb_min_child_weight" = 1),
                                     nthread = 8, # number of cores, you can change,
                                     covar_bl_method = "absolute",
                                     covar_bl_trs = 0.1,
                                     trim_quantiles = c(0.01,0.99), # trimed, you can change,
                                     optimized_compile = TRUE, #created a column counter for how many times matched,
                                     max_attempt = 1,
                                     matching_fun = "matching_l1",
                                     delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                     scale = 1.0)
# Best Mean absolute correlation:  0.0524237671382134 | Covariate balance threshold:  0.1
match_pop_data <- match_pop_all$pseudo_pop
colnames(match_pop_data)[1] <- "zip"

aggregate_data_rm_counter <- merge(aggregate_data_rm, match_pop_data[, c("year", "zip", "counter")], by = c("year", "zip"), all.y = TRUE)
aggregate_data_rm_counter <- subset(aggregate_data_rm_counter, counter > 0)

matchingrm_gnm <- summary(gnm(dead ~ pm25 + offset(log(time_count)), 
                              eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
                              data = aggregate_data_rm_counter,
                              family = poisson(link = "log"),
                              weights = counter))
exp(10*matchingrm_gnm$coefficients[1])
# [1] 1.036226

# trimed
# GPS: PM | potential confounders + ozone
match_pop_ozone <- generate_pseudo_pop(Y = zip_year$zip,
                                     w = zip_year$pm25,
                                     c = zip_year[, c(2, 10:13,15:16,18)],
                                     ci_appr = "matching",
                                     pred_model = "sl",
                                     gps_model = "parametric",
                                     use_cov_transform = FALSE,
                                     #transformers = list("pow2", "pow3"),
                                     sl_lib = c("m_xgboost"),
                                     params = list("xgb_nrounds" = 50,
                                                   "xgb_max_depth" = 6,
                                                   "xgb_eta" = 0.3,
                                                   "xgb_min_child_weight" = 1),
                                     nthread = 8, # number of cores, you can change,
                                     covar_bl_method = "absolute",
                                     covar_bl_trs = 0.1,
                                     trim_quantiles = c(0.01,0.99), # trimed, you can change,
                                     optimized_compile = TRUE, #created a column counter for how many times matched,
                                     max_attempt = 1,
                                     matching_fun = "matching_l1",
                                     delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                     scale = 1.0)
# Best Mean absolute correlation:  0.0456610146068809 | Covariate balance threshold:  0.1
match_pop_ozone_data <- match_pop_ozone$pseudo_pop
colnames(match_pop_ozone_data)[1] <- "zip"

aggregate_data_rm_ozone_counter <- merge(aggregate_data_rm, match_pop_ozone_data[, c("year", "zip", "counter")], by = c("year", "zip"), all.y = TRUE)
aggregate_data_rm_ozone_counter <- subset(aggregate_data_rm_ozone_counter, counter > 0)

matchingrm_gnm_ozone <- summary(gnm(dead ~ pm25 + offset(log(time_count)), 
                              eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
                              data = aggregate_data_rm_ozone_counter,
                              family = poisson(link = "log"),
                              weights = counter))
exp(10*matchingrm_gnm_ozone$coefficients[1])

# trimed
# GPS: PM | potential confounders + no2
#c2 = as.data.frame(data.matrix(zip_year[, c(2, 10:13,15:16,19)]))
match_pop_no2 <- generate_pseudo_pop(Y = zip_year$zip,
                                       w = zip_year$pm25,
                                       c = zip_year[, c(2, 10:13,15:16,19)],
                                       ci_appr = "matching",
                                       pred_model = "sl",
                                       gps_model = "parametric",
                                       use_cov_transform = FALSE,
                                       #transformers = list("pow2", "pow3"),
                                       sl_lib = c("m_xgboost"),
                                       params = list("xgb_nrounds" = 50,
                                                     "xgb_max_depth" = 6,
                                                     "xgb_eta" = 0.3,
                                                     "xgb_min_child_weight" = 1),
                                       nthread = 8, # number of cores, you can change,
                                       covar_bl_method = "absolute",
                                       covar_bl_trs = 0.1,
                                       trim_quantiles = c(0.01,0.99), # trimed, you can change,
                                       optimized_compile = TRUE, #created a column counter for how many times matched,
                                       max_attempt = 1,
                                       matching_fun = "matching_l1",
                                       delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                       scale = 1.0)
# Best Mean absolute correlation:  0.0581104418500892 | Covariate balance threshold:  0.1
match_pop_no2_data <- match_pop_no2$pseudo_pop
colnames(match_pop_no2_data)[1] <- "zip"
#match_pop_no2_data$year = as.factor(levels(zip_year$year)[match_pop_no2_data$year])
#match_pop_no2_data$zip = as.factor(match_pop_no2_data$zip)
aggregate_data_rm_no2_counter <- merge(aggregate_data_rm, match_pop_no2_data[, c("year", "zip", "counter")], by = c("year", "zip"), all.y = TRUE)
aggregate_data_rm_no2_counter <- subset(aggregate_data_rm_no2_counter, counter > 0)

matchingrm_gnm_no2_pilot <- gnm(dead ~ pm25,
                                  offset = log(time_count), 
                                  eliminate = as.factor(sex):as.factor(dual):as.factor(race):as.factor(entry_age_break):as.factor(followup_year),
                                  data = subset(aggregate_data_rm_no2_counter, entry_age_break != 8),
                                  family = poisson(link = "log"),
                                  weights = counter,
                                  trace = TRUE,
                                  tolerance = 1e-10)
exp(10*matchingrm_gnm_no2_pilot$coefficients[1])

cont_table_all = table(aggregate_data_rm_no2_counter2$sex, 
                    aggregate_data_rm_no2_counter2$dual, 
                    aggregate_data_rm_no2_counter2$race, 
                    aggregate_data_rm_no2_counter2$entry_age_break, 
                    aggregate_data_rm_no2_counter2$followup_year)

aggregate_data_rm_no2_counter2 = aggregate_data_rm_no2_counter
#aggregate_data_rm_no2_counter2$followup_year[aggregate_data_rm_no2_counter2$followup_year >= 12] = 12
#aggregate_data_rm_no2_counter2$entry_age_break[aggregate_data_rm_no2_counter2$entry_age_break >= 7] = 7
#aggregate_data_rm_no2_counter2$race[aggregate_data_rm_no2_counter2$race == 0] = 3
aggregate_data_rm_no2_counter2 = subset(aggregate_data_rm_no2_counter2, 
                                        race != 0)

matchingrm_gnm_no2 <- summary(gnm(dead ~ pm25,
                                    offset = log(time_count), 
                                    eliminate = as.factor(sex):as.factor(dual):as.factor(race):as.factor(entry_age_break):as.factor(followup_year),
                                    data = aggregate_data_rm_no2_counter2,
                                    family = poisson(link = "log"),
                                    weights = counter,
                                  trace = TRUE))
exp(10*matchingrm_gnm_no2$coefficients[1])
# 0.9909953

# trimed
# GPS: PM | potential confounders + ox
match_pop_ox <- generate_pseudo_pop(Y = zip_year$zip,
                                     w = zip_year$pm25,
                                     c = zip_year[, c(2, 10:13,15:16,20)],
                                     ci_appr = "matching",
                                     pred_model = "sl",
                                     gps_model = "parametric",
                                     use_cov_transform = FALSE,
                                     #transformers = list("pow2", "pow3"),
                                     sl_lib = c("m_xgboost"),
                                     params = list("xgb_nrounds" = 50,
                                                   "xgb_max_depth" = 6,
                                                   "xgb_eta" = 0.3,
                                                   "xgb_min_child_weight" = 1),
                                     nthread = 8, # number of cores, you can change,
                                     covar_bl_method = "absolute",
                                     covar_bl_trs = 0.1,
                                     trim_quantiles = c(0.01,0.99), # trimed, you can change,
                                     optimized_compile = TRUE, #created a column counter for how many times matched,
                                     max_attempt = 1,
                                     matching_fun = "matching_l1",
                                     delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                     scale = 1.0)
match_pop_ox_data <- match_pop_ox$pseudo_pop
colnames(match_pop_ox_data)[1] <- "zip"

aggregate_data_rm_ox_counter <- merge(aggregate_data_rm, match_pop_ox_data[, c("year", "zip", "counter")], by = c("year", "zip"), all.y = TRUE)
aggregate_data_rm_ox_counter <- subset(aggregate_data_rm_ox_counter, counter > 0)

aggregate_data_rm_ox_counter2 = aggregate_data_rm_ox_counter
aggregate_data_rm_ox_counter2$followup_year[aggregate_data_rm_ox_counter2$followup_year >= 12] = 12
aggregate_data_rm_ox_counter2$entry_age_break[aggregate_data_rm_ox_counter2$entry_age_break >= 7] = 7
aggregate_data_rm_ox_counter2 = subset(aggregate_data_rm_ox_counter2, !(race %in% c(0,6)))

cont_table_all = table(aggregate_data_rm_ox_counter2$sex, 
                       aggregate_data_rm_ox_counter2$dual, 
                       aggregate_data_rm_ox_counter2$race, 
                       aggregate_data_rm_ox_counter2$entry_age_break, 
                       aggregate_data_rm_ox_counter2$followup_year)

matchingrm_gnm_ox <- summary(gnm(dead ~ pm25 + ox + offset(log(time_count)), 
                                  eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
                                  data = aggregate_data_rm_ox_counter2,
                                  family = poisson(link = "log"),
                                  weights = counter,
                                 trace = TRUE,
                                 tolerance = 1e-10))
exp(10*matchingrm_gnm_ox$coefficients[1])
## no trim
## Error in if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) { : 
## missing value where TRUE/FALSE needed
match_pop_all <- generate_pseudo_pop(Y = zip_year$zip,
                                     w = zip_year$pm25,
                                     c = zip_year[, c(2, 10:13,15:16)],
                                     ci_appr = "matching",
                                     pred_model = "sl",
                                     gps_model = "parametric",
                                     use_cov_transform = FALSE,
                                     #transformers = list("pow2", "pow3"),
                                     sl_lib = c("m_xgboost"),
                                     params = list("xgb_nrounds" = 50,
                                                   "xgb_max_depth" = 6,
                                                   "xgb_eta" = 0.3,
                                                   "xgb_min_child_weight" = 1),
                                     nthread = 8, # number of cores, you can change,
                                     covar_bl_method = "absolute",
                                     covar_bl_trs = 0.1,
                                     trim_quantiles = c(0.00,1.00), # trimed, you can change,
                                     optimized_compile = TRUE, #created a column counter for how many times matched,
                                     max_attempt = 1,
                                     matching_fun = "matching_l1",
                                     delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                     scale = 1.0)
match_pop_data <- match_pop_all$pseudo_pop
colnames(match_pop_data)[1] <- "zip"

aggregate_data_rm_counter <- merge(aggregate_data_rm, match_pop_data[, c("year", "zip", "counter")], by = c("year", "zip"), all.y = TRUE)
aggregate_data_rm_counter <- subset(aggregate_data_rm_counter, counter > 0)

matchingrm_gnm <- summary(gnm(dead ~ pm25 + offset(log(time_count)), 
                              eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
                              data = aggregate_data_rm_counter,
                              family = poisson(link = "log"),
                              weights = counter,
                              trace = TRUE,
                              tolerance = 1e-20))
exp(10*matchingrm_gnm$coefficients[1])



