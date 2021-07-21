* Written by R;
*  write.foreign(national_merged2016, datafile = "/nfs/home/X/xwu/shared_space/ci3_health_data/medicare/mortality/1999_2016/wu/output_data/harmonized_national2016.csv",  ;

DATA  national2016 ;
LENGTH
 zip $ 5
 dual $ 1
 statecode $ 2
 region $ 9
;

INFILE  "/nfs/home/X/xwu/shared_space/ci3_health_data/medicare/mortality/1999_2016/wu/output_data/harmonized_national2016.csv" 
     DSD 
     LRECL= 184 ;
INPUT
 year
 zip
 sex
 race
 age
 dual
 entry_age_break
 entry_age
 statecode
 followup_year
 followup_year_plus_one
 dead
 medhouseholdincome
 medianhousevalue
 poverty
 education
 popdensity
 pct_owner_occ
 region
 pm25
;
RUN;


ODS HTML FILE='model1_2016.html' path='/nfs/home/X/xwu/shared_space/ci3_analysis/HEI_Final/Harmonized/';PROC phreg data = national2016;
class sex year / ref=first;
model (followup_year,followup_year_plus_one)* dead(0) = 
pm25 entry_age sex year/rl ties=efron;
title "model1; 2000-2016";
run;ODS HTML CLOSE;
