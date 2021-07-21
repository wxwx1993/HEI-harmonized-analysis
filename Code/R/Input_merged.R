library(fst)
library(data.table)
library(parallel)
library("dplyr")
library("foreign")
library("gnm")

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

# write .csv file for SAS Cox analysis
write.foreign(national_merged2016, datafile="/nfs/home/X/xwu/shared_space/ci3_health_data/medicare/mortality/1999_2016/wu/output_data/harmonized_national2016.csv", 
              codefile="/nfs/home/X/xwu/shared_space/ci3_analysis/HEI_Final/Harmonized/harmonized_national2016.sas",
              package = c("SAS"))

# aggregate data for R Poisson analysis
national_merged2016$time_count<-national_merged2016$followup_year_plus_one-national_merged2016$followup_year
dead_personyear<-aggregate(cbind(national_merged2016$dead, national_merged2016$time_count),
                           by=list(national_merged2016$zip,
                                   national_merged2016$year,
                                   national_merged2016$sex,
                                   national_merged2016$race,
                                   national_merged2016$dual,
                                   national_merged2016$entry_age_break,
                                   national_merged2016$followup_year), FUN=sum)

confounders<-aggregate(national_merged2016[,c(13:20)], 
                       by=list(national_merged2016$zip,
                               national_merged2016$year,
                               national_merged2016$sex,
                               national_merged2016$race,
                               national_merged2016$dual,
                               national_merged2016$entry_age_break,
                               national_merged2016$followup_year), 
                       FUN=min)

aggregate_data_rm<-merge(dead_personyear,confounders
                         ,by=c("Group.1", "Group.2", "Group.3", "Group.4", "Group.5", "Group.6", "Group.7"))
colnames(aggregate_data_rm)[8:9]<-c("dead","time_count")
colnames(aggregate_data_rm)[1:7]<-c("zip","year","sex","race","dual","entry_age_break","followup_year")
aggregate_data_rm<-subset(aggregate_data_rm[complete.cases(aggregate_data_rm) ,])

save(aggregate_data_rm, file="/nfs/nsaph_ci3/ci3_analysis/HEI_Final/Harmonized/aggregate_data_rm.RData")
