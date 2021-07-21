

ODS HTML FILE='model3_2016.html' path='/nfs/home/X/xwu/shared_space/ci3_analysis/HEI_Final/Harmonized/';PROC phreg data = national2016;
class sex year dual race/ ref=first;
model (followup_year,followup_year_plus_one)* dead(0) = 
pm25 entry_age sex year dual race 
medhouseholdincome medianhousevalue poverty education pct_owner_occ/rl ties=efron;
title "model3; 2000-2016";
run;ODS HTML CLOSE;
