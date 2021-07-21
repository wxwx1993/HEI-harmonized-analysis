# Bootstrap on zip-code cluster to obtain robust CIs account for spatial correlation
library("mgcv")
library("survival")
library("parallel")
library("dplyr")
library(fst)
library(data.table)

f <- list.files("/nfs/nsaph_ci3/ci3_health_data/medicare/mortality/1999_2016/wu/cache_data/merged_by_year_v2",
                pattern = "\\.fst",
                full.names = TRUE)

myvars <- c("year","zip","sex","race","age","dual","entry_age_break","entry_age","statecode",
            "followup_year","followup_year_plus_one","dead",
            "medhouseholdincome","medianhousevalue",
            "poverty","education","popdensity", "pct_owner_occ")
national_merged2016 <- rbindlist(lapply(f,
                                        read_fst,
                                        columns = myvars,
                                        as.data.table=FALSE))

national_merged2016$zip <- sprintf("%05d", national_merged2016$zip)

NORTHEAST=c("NY","MA","PA","RI","NH","ME","VT","CT","NJ")  
SOUTH=c("DC","VA","NC","WV","KY","SC","GA","FL","AL","TN","MS","AR","MD","DE","OK","TX","LA")
MIDWEST=c("OH","IN","MI","IA","MO","WI","MN","SD","ND","IL","KS","NE")
WEST=c("MT","CO","WY","ID","UT","NV","CA","OR","WA","AZ","NM")

national_merged2016$region=ifelse(national_merged2016$state %in% NORTHEAST, "NORTHEAST", 
                                  ifelse(national_merged2016$state %in% SOUTH, "SOUTH",
                                         ifelse(national_merged2016$state %in% MIDWEST, "MIDWEST",
                                                ifelse(national_merged2016$state %in% WEST, "WEST",
                                                       NA))))
pm_rm<-NA
for(i in 2000:2016){
  temp<-read.csv(paste0("/nfs/home/X/xwu/shared_space/ci3_exposure/pm25/whole_us/annual/zipcode/rm_predictions/ben_2019_10_29/data_acag_pm25_zip-year/zip_pm25_",i, ".csv"))
  temp$YEAR<-rep(i, nrow(temp))
  pm_rm<-rbind(pm_rm, temp)
  
}
rm(temp)
pm_rm$ZIP <- sprintf("%05d", pm_rm$ZIP)

national_merged2016 <- merge(national_merged2016, 
                             pm_rm[, c(1,3:4)], 
                             by.x = c("year", "zip"),
                             by.y = c("YEAR", "ZIP"),
                             all.x = TRUE)

national_merged2016 <- national_merged2016[complete.cases(national_merged2016), ]


load("/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/aggregate_data_rm.RData")
all_zip <- unique(aggregate_data_rm$zip)
rm(aggregate_data_rm)
num_uniq_zip <- length(all_zip)

# Save the bootstrapped data to accelarate computing
dir.create(file.path("/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/Boots/Cox_boots"), showWarnings = FALSE)
lapply(455:500,function(boots_id){
  set.seed(boots_id)
  zip_sample<-sample(1:num_uniq_zip,floor(2*sqrt(num_uniq_zip)),replace=T) 
  national_merged2016_boots <- subset(national_merged2016, zip %in% all_zip[zip_sample]) 
  write_fst(national_merged2016_boots, paste0("/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/Boots/Cox_boots/",
                        boots_id,".fst"))
})

## mod1
Cox_coefs_boots<-NULL
Cox_coefs_boots <- sapply(1:500, function(boots_id) {
  set.seed(boots_id)
  national_merged2016_boots<- read_fst(paste0("/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/Boots/Cox_boots/",
                                       boots_id,".fst"))
    ######## Main models ########
  Cox <- coxph(Surv(followup_year, followup_year_plus_one, dead) ~ pm25 +
                 (entry_age) + as.factor(sex) + as.factor(year),
               data=national_merged2016_boots,
               ties = c("efron"),
               na.action = na.omit)
  return(summary(Cox)$coefficients[1])
})

save(num_uniq_zip,
     Cox_coefs_boots,
     file="/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/Boots/Cox_coefs_boots_mod1.RData")

exp(log(1.003)+10*1.96*sd(Cox_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip))
exp(log(1.003)-10*1.96*sd(Cox_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip))
#> exp(mean(Cox_coefs_boots)*10)
#[1] 1.005475

## mod2
Cox_coefs_boots<-NULL
Cox_coefs_boots <- sapply(1:500, function(boots_id) {
  set.seed(boots_id)
  national_merged2016_boots<- read_fst(paste0("/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/Boots/Cox_boots/",
                                              boots_id,".fst"))
  ######## Main models ########
  Cox <- coxph(Surv(followup_year, followup_year_plus_one, dead) ~ pm25 +
                 (entry_age) + as.factor(sex) + as.factor(year) +
                 + as.factor(dual) + as.factor(race),
               data=national_merged2016_boots,
               ties = c("efron"),
               na.action = na.omit)
  return(summary(Cox)$coefficients[1])
})

save(num_uniq_zip,
     Cox_coefs_boots,
     file="/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/Boots/Cox_coefs_boots_mod2.RData")

exp(log(1.013)+10*1.96*sd(Cox_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip))
exp(log(1.013)-10*1.96*sd(Cox_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip))

## mod3
Cox_coefs_boots<-NULL
Cox_coefs_boots <- sapply(1:500, function(boots_id) {
  set.seed(boots_id)
  national_merged2016_boots<- read_fst(paste0("/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/Boots/Cox_boots/",
                                              boots_id,".fst"))
  ######## Main models ########
  Cox <- coxph(Surv(followup_year, followup_year_plus_one, dead) ~ pm25 +
                 (entry_age) + as.factor(sex) + as.factor(year) +
                 + as.factor(dual) + as.factor(race) +
                 medhouseholdincome + medianhousevalue + poverty + education + pct_owner_occ,
               data=national_merged2016_boots,
               ties = c("efron"),
               na.action = na.omit)
  return(summary(Cox)$coefficients[1])
})

save(num_uniq_zip,
     Cox_coefs_boots,
     file="/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/Boots/Cox_coefs_boots_mod3.RData")

exp(log(1.054)+10*1.96*sd(Cox_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip))
exp(log(1.054)-10*1.96*sd(Cox_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip))

## mod4
Cox_coefs_boots<-NULL
Cox_coefs_boots <- sapply(1:500, function(boots_id) {
  set.seed(boots_id)
  national_merged2016_boots<- read_fst(paste0("/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/Boots/Cox_boots/",
                                              boots_id,".fst"))
  ######## Main models ########
  Cox <- coxph(Surv(followup_year, followup_year_plus_one, dead) ~ pm25 +
                 (entry_age) + as.factor(sex) + as.factor(year) +
                 as.factor(dual) + as.factor(race) + 
                 medhouseholdincome + medianhousevalue + poverty + education + pct_owner_occ +
                 as.factor(region),
               data=national_merged2016_boots,
               ties = c("efron"),
               na.action = na.omit)
  return(summary(Cox)$coefficients[1])
})

save(num_uniq_zip,
     Cox_coefs_boots,
     file="/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/Boots/Cox_coefs_boots_mod4.RData")

exp(log(1.051)+10*1.96*sd(Cox_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip))
exp(log(1.051)-10*1.96*sd(Cox_coefs_boots) *sqrt(2*sqrt(num_uniq_zip))/sqrt(num_uniq_zip))

