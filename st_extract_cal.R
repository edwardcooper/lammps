
source('~/Dropbox/lammps/equilibrium_extract.R', echo=TRUE)
 
#library(beepr) # for sound effect.
 log_data=read.table("log320_solomon.csv")
 log_data=log_data%>%long_to_wide()
 colnames(log_data)
 library(purrr)
 library(magrittr)
 st_test_result=log_data%>%map(function(data) st_extract_tc(data,method="kpss"))
 st_test_result

 
 # library(purrr)
 # library(magrittr)
 # st_test_result=log_data%>%map(function(data) st_extract_tc(data,method="kpss"))
 # st_test_result
 
 
log_data3=read.table("log.csv")
log_data3=log_data3%>%long_to_wide()
colnames(log_data3)
st_test_result3=log_data3%>%map(function(data) st_extract_tc(data,method="fractal"))
st_test_result4=log_data3%>%map(function(data) st_extract_tc(data,method="adf"))
st_test_result5=log_data3%>%map(function(data) st_extract_tc(data,method="kpss"))
st_test_result6=log_data3%>%map(function(data) st_extract_tc(data,method="pp"))
st_test_result7=log_data3%>%map(function(data) st_extract_tc(data,method="stats.pp"))

st_test_result3
st_test_result4
st_test_result5
st_test_result6
st_test_result7