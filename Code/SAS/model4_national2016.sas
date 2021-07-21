

ODS HTML FILE='model4_2016.html' path='/nfs/home/X/xwu/shared_space/ci3_analysis/HEI_Final/Harmonized/';PROC phreg data = national2016;
class sex year dual race region/ ref=first;
model (followup_year,followup_year_plus_one)* dead(0) = 
pm25 entry_age sex year dual race 
medhouseholdincome medianhousevalue poverty education pct_owner_occ
region/rl ties=efron;
title "model4; 2000-2016";
run;ODS HTML CLOSE;
