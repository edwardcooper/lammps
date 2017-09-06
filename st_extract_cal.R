
source('~/Dropbox/lammps/equilibrium_extract.R', echo=TRUE)
 
#library(beepr) # for sound effect.
 log_data=read.table("log320_solomon.csv")
 log_data=log_data%>%long_to_wide()
 colnames(log_data)
 library(purrr)
 library(magrittr)
 st_test_result=log_data%>%map(function(data) st_extract_tc(data,method="kpss"))
 st_test_result%>%map(length)

 
 # library(purrr)
 # library(magrittr)
 # st_test_result=log_data%>%map(function(data) st_extract_tc(data,method="kpss"))
 # st_test_result
log_data4=read.table("log600K.csv")
log_data4=log_data4%>%long_to_wide()
st_test_result8=log_data4%>%map(function(data) st_extract_tc(data,method="fractal",interval=10^5))
st_test_result9=log_data4%>%map(function(data) st_extract_tc(data,method="kpss"))

st_test_result8%>%map(length)
st_test_result9%>%map(length)

log_data3=read.table("log.csv")
log_data3=log_data3%>%long_to_wide()
colnames(log_data3)
st_test_result3=log_data3%>%map(function(data) st_extract_tc(data,method="fractal",interval=10^6))
#st_test_result4=log_data3%>%map(function(data) st_extract_tc(data,method="adf"))
st_test_result5=log_data3%>%map(function(data) st_extract_tc(data,method="kpss"))
#st_test_result6=log_data3%>%map(function(data) st_extract_tc(data,method="pp"))
#st_test_result7=log_data3%>%map(function(data) st_extract_tc(data,method="stats.pp",interval=10^6))

st_test_result3%>%map(length)
st_test_result4%>%map(length)
st_test_result5%>%map(length)
st_test_result6%>%map(length)
st_test_result7%>%map(length)
