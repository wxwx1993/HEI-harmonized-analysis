# Bootstrap on zip-code cluster to obtain robust CIs account for spatial correlation
library("mgcv")
library("gnm")
library("parallel")
library("dplyr")
require(parallel)
library(data.table)


dir_data = '/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/'
dir_out = '/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/Boots/'

load(paste0(dir_data,"aggregate_data_rm.RData"))

# poisson
aggregate_data.list<-split(aggregate_data_rm, list(aggregate_data_rm$zip))
num_uniq_zip <- length(unique(aggregate_data_rm$zip))

#mod1
loglinear_coefs_boots <- mclapply(1:500, function(boots_id) {
  set.seed(boots_id)
  zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)), replace = TRUE) 
  aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
  
  gnm_pm25_mod1 <- gnm(dead ~ pm25 + as.factor(year) +
                         offset(log(time_count)),
                       eliminate = (as.factor(sex):as.factor(entry_age_break):as.factor(followup_year)),
                       data = aggregate_data_boots,
                       family = poisson(link = "log"))
  
  return(gnm_pm25_mod1$coefficients[1])
  rm(aggregate_data_boots)
}, mc.cores = 8)

save(num_uniq_zip,
     loglinear_coefs_boots,
     file="/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/Boots/loglinear_coefs_boots_mod1.RData")

exp(10*(summary_pm25_mod1$coefficients[1]-1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))
exp(10*(summary_pm25_mod1$coefficients[1]+1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))

#mod2
loglinear_coefs_boots <- sapply(1:500, function(boots_id) {
  set.seed(boots_id)
  zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)), replace = TRUE) 
  aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
  
  gnm_pm25_mod2 <- gnm(dead ~ pm25 + as.factor(year) +
                         offset(log(time_count)),
                       eliminate = (as.factor(sex):as.factor(entry_age_break):as.factor(race):as.factor(dual):as.factor(followup_year)),
                       data = aggregate_data_boots,
                       family = poisson(link = "log"))
  
  return(gnm_pm25_mod2$coefficients[1])
  rm(aggregate_data_boots)
})

save(num_uniq_zip,
     loglinear_coefs_boots,
     file="/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/Boots/loglinear_coefs_boots_mod2.RData")

exp(10*(summary_pm25_mod2$coefficients[1]-1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))
exp(10*(summary_pm25_mod2$coefficients[1]+1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))

#mod3
loglinear_coefs_boots <- sapply(1:500, function(boots_id) {
  set.seed(boots_id)
  zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)), replace = TRUE) 
  aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
  
  gnm_pm25_mod3 <- gnm(dead ~ pm25 + as.factor(year) + 
                         medhouseholdincome + medianhousevalue + poverty + education + pct_owner_occ +
                         offset(log(time_count)),
                       eliminate = (as.factor(sex):as.factor(entry_age_break):as.factor(race):as.factor(dual):as.factor(followup_year)),
                       data = aggregate_data_boots,
                       family = poisson(link = "log"))
  
  return(gnm_pm25_mod3$coefficients[1])
  rm(aggregate_data_boots)
})

save(num_uniq_zip,
     loglinear_coefs_boots,
     file="/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/Boots/loglinear_coefs_boots_mod3.RData")

exp(10*(summary_pm25_mod3$coefficients[1]-1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))
exp(10*(summary_pm25_mod3$coefficients[1]+1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))

#mod4
loglinear_coefs_boots <- mclapply(1:500, function(boots_id) {
  set.seed(boots_id)
  zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)), replace = TRUE) 
  aggregate_data_boots<-data.frame(Reduce(rbind,aggregate_data.list[zip_sample]))
  
  gnm_pm25_mod4 <- gnm(dead ~ pm25 + as.factor(year) +
                         medhouseholdincome + medianhousevalue + poverty + education + pct_owner_occ +
                         as.factor(region) + 
                         offset(log(time_count)),
                       eliminate = (as.factor(sex):as.factor(entry_age_break):as.factor(race):as.factor(dual):as.factor(followup_year)),
                       data = aggregate_data_boots,
                       family = poisson(link = "log"))
  return(gnm_pm25_mod4$coefficients[1])
  rm(aggregate_data_boots)
},mc.cores = 8)

save(num_uniq_zip,
     loglinear_coefs_boots,
     file="/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/Boots/loglinear_coefs_boots_mod4.RData")

exp(10*(summary_pm25_mod4$coefficients[1]-1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))
exp(10*(summary_pm25_mod4$coefficients[1]+1.96*sd(loglinear_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip)))
