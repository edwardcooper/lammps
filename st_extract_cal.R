
source("https://raw.githubusercontent.com/edwardcooper/lammps/master/equilibrium_extract.R", echo=TRUE)
source("https://raw.githubusercontent.com/edwardcooper/lammps/master/map_pc.R", echo=TRUE) 
source("https://raw.githubusercontent.com/edwardcooper/mlmodel_select/master/timeRecord_functions.R",echo=TRUE)


#library(beepr) # for sound effect.
timeRecordB()
library(purrr)
library(magrittr)
log_data320k=read.table("log320K.csv")
timeRecordB(output_message="Load_data")
log_data320k=log_data320k%>%long_to_wide()
timeRecordB(output_message="long_to_wide")
colnames(log_data320k)
timeRecordB()
st_test_result320k=log_data320k%>%map_pc(function(data) st_extract_tc(data,method="kpss",interval=nrow(log_data320k)/10))
st_test_result320k%>%map(length)

timeRecordB(output_message="test_calculation")

timeRecordR()



 library(purrr)
 library(magrittr)
 timeRecordB()
 log_data600k=read.table("log600K.csv")
 timeRecordB(output_message="Load_data")
 
 log_data600k=log_data600k%>%long_to_wide()
 timeRecordB(output_message="long_to_wide")
 colnames(log_data600k)
 timeRecordB()
 st_test_result600k=log_data600k%>%map_pc(function(data) st_extract_tc(data,method="kpss",interval=nrow(log_data320k)/10) )
 st_test_result600k%>%map(length)
 timeRecordB(output_message="test_calculation")
 
 
timeRecordR()%>%filter(output_message!="None")%>%select(output_message,run_time) 


#######################################################################################################################

timeRecordB()
st_test_result320k=log_data320k%>%map_pc(function(data) st_extract_tc(data,method="kpss",interval=nrow(log_data320k)/10),nthread=1,type="FORK")
st_test_result320k%>%map(length)

timeRecordB(output_message="test_calculation core_num: 1")
#######################################################################################################################
timeRecordB()
st_test_result320k=log_data320k%>%map_pc(function(data) st_extract_tc(data,method="kpss",interval=nrow(log_data320k)/10),nthread=2,type="FORK")
st_test_result320k%>%map(length)
timeRecordB(output_message="test_calculation core_num: 2")
#######################################################################################################################
timeRecordB()
st_test_result320k=log_data320k%>%map_pc(function(data) st_extract_tc(data,method="kpss",interval=nrow(log_data320k)/10),nthread=3,type="FORK")
st_test_result320k%>%map(length)

timeRecordB(output_message="test_calculation core_num: 3")
#######################################################################################################################

timeRecordB()
st_test_result320k=log_data320k%>%map_pc(function(data) st_extract_tc(data,method="kpss",interval=nrow(log_data320k)/10),nthread=4,type="FORK")
st_test_result320k%>%map(length)

timeRecordB(output_message="test_calculation core_num: 4")
#######################################################################################################################
timeRecordB()
st_test_result320k=log_data320k%>%map_pc(function(data) st_extract_tc(data,method="kpss",interval=nrow(log_data320k)/10),nthread=5,type="FORK")
st_test_result320k%>%map(length)

timeRecordB(output_message="test_calculation core_num: 5")
#######################################################################################################################


timeRecordB()
st_test_result320k=log_data320k%>%map_pc(function(data) st_extract_tc(data,method="kpss",interval=nrow(log_data320k)/10),nthread=6,type="FORK")
st_test_result320k%>%map(length)

timeRecordB(output_message="test_calculation core_num: 6")
#######################################################################################################################

timeRecordB()
st_test_result320k=log_data320k%>%map_pc(function(data) st_extract_tc(data,method="kpss",interval=nrow(log_data320k)/10),nthread=7,type="FORK")
st_test_result320k%>%map(length)

timeRecordB(output_message="test_calculation core_num: 7")
#######################################################################################################################

timeRecordB()
st_test_result320k=log_data320k%>%map_pc(function(data) st_extract_tc(data,method="kpss",interval=nrow(log_data320k)/10),nthread=8,type="FORK")
st_test_result320k%>%map(length)

timeRecordB(output_message="test_calculation core_num: 8")
#######################################################################################################################





timeRecordR()%>%filter(output_message!="None")%>%select(output_message,run_time) 


