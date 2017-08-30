# The purpose of this function is to get the part of simulation that is in equilibrium from the energies, temperatures and etc. 

library(magrittr)



long_to_wide=function(dataset){
  
  # change the colnames of dataset for easier manipulation.
  colnames(dataset)=c("property","value")
  
  # Find out all the available values in the dataset.
  dataset_avail=table(dataset$property)%>%as.data.frame()
  
  # extract the names and convert it into characters vector.
  names=dataset_avail$Var1%>%as.character()
  
  # I use foreach here to allow parallel processing if necessary in the future.
  library(foreach)
  dataset=foreach(i=1:length(names),.combine = cbind, .packages = "dplyr" )%do%{
    
    #column_value=dataset%>%filter(property==names[i])
    column_value=dataset%>%filter(property==names[i])%>%select(value)
    return(column_value)
  }
  
  colnames(dataset)=names
  #dataset=cbind(timestep=1:dim(dataset)[1],dataset)
  gc()
  return(dataset)
  
}


#################################################
# main function in this script




# Could only chosee alpha in (0.01, 0.1)
st_extract_spectra=function(data, del_num=1000,alpha=0.05){
  # do the statistical test on the data
  test_fractal=fractal::stationarity(data) #null hypothesis: stationary. 
  # extract the p value. 
  p_value=attr(test_fractal,"pvals")[1]
  # test to see if the p value is smaller than alpha (type I error rate). 
  # If it is TRUE, it means that the p value is smaller than alpha, then null hypothesis that the data is staionary is rejected. Data not stationary. 
  # If it is FALSE, it means that the p value is bigger than alpha, then we fail to reject the null hypothesis that the data is staionary. Data stationary. 
  
  
  non_stationary=p_value<alpha
  
  # If it is stationary, then we will return the original data.
  # If it is not stationary, then we will need to chop off some of the data to test again until we find the part of data that is stationary. 
  if(!non_stationary){
    
    result_data=data
    plot(result_data,main="Entire stationary data")
    print("This data is stationary")
  }else if (non_stationary){
    
    # while non_stationary is 
    while(non_stationary){
      del_index=seq(from=1,to=del_num,by=1)
      #print(length(del_index))
      data<-data[-del_index]
      test_fractal=fractal::stationarity(data)
      p_value=attr(test_fractal,"pvals")[1]
      print(p_value)
      non_stationary=p_value<alpha
      #print(length(data))
    }
    
    result_data=data
    plot(result_data,main="stationary part")
    #print("The length of the data is", length(data))
  }
  gc()
  
  return(result_data)
}

# test use
#rnorm(100)%>%st_extract_spectra()
###############################################################################

library(beepr)
log_data=read.table("log_npt.csv")
log_data=log_data%>%long_to_wide()
colnames(log_data)

###############################################################################

# Add error handling, it will return NA value if the data has no staionary part with the precision defined.
st_extract_spectra_tc=function(data){
  out=tryCatch(st_extract_spectra(data,alpha=0.05,del_num = 10000), error=function(e) { message("The data has no stationary part with the precision defined.") return(NA) } )
  return(out)
}

# test usage for the function st_extract_spectra_tc
#KinEng_eq=st_extract_spectra_tc(log_data[,3]) # not stationary

# Do the st_extract_spectra_tc on all the features. 
st_test_result=apply(log_data,2,st_extract_spectra_tc)
###############################################################################


# KinEng_eq=log_data[,2]%>%st_extract_spectra_tc(alpha=0.05,del_num = 10000) # not stationary
# 
# PotEng_eq=log_data[,3]%>%st_extract_spectra(alpha=0.05,del_num = 10000) # not stationary
# 
# 
# Press_eq=log_data[,4]%>%st_extract_spectra(alpha=0.05,del_num = 10000) # not stationary
# print("This is pressure data.")
# beep(8)
# Temp_eq=log_data[,5]%>%st_extract_spectra(alpha=0.05,del_num = 10000) # not stationary
# print("This is temperature data.")
# beep(8)
# TotEng_eq=log_data[,6]%>%st_extract_spectra(alpha=0.05,del_num = 10000) # Staionary for the last 10001 data points
# print("This is total energy data")
# beep(8)
# Volume_eq=log_data[,7]%>%st_extract_spectra(alpha=0.05,del_num = 10000) 
# print("This is volume data")
# beep(8)
# # test_kinEng=fractal::stationarity(KinEng_eq)
# # attr(test_kinEng,"pvals")[1]>0.05

################################################################################


################################################################################
# library(beepr)
# log_solomon320=read.table("log320_solomon.csv")
# log_solomon320=log_solomon320%>%long_to_wide()
# colnames(log_solomon320)
# 
# 
# KinEng_eq=log_solomon320[,2]%>%st_extract_spectra(alpha=0.05,del_num = 50) # stationary 
# print("This is Kinetic energy data")
# 
# PotEng_eq=log_solomon320[,3]%>%st_extract_spectra(alpha=0.05,del_num = 50) # stationary 
# print("This is Potential energy data")
# Press_eq=log_solomon320[,4]%>%st_extract_spectra(alpha=0.05,del_num = 50)# stationary 
# print("This is pressure data.")
# beep(8)
# Temp_eq=log_solomon320[,5]%>%st_extract_spectra(alpha=0.05,del_num = 50) # stationary 
# print("This is temperature data.")
# beep(8)
# TotEng_eq=log_solomon320[,6]%>%st_extract_spectra(alpha=0.05,del_num = 50) # stationary after 50 data points 
# print("This is total energy data")
# beep(8)
#######################################################################################

#######################################################################################
# library(beepr)
# log_solomon345=read.table("log345_solomon.csv")
# log_solomon345=log_solomon345%>%long_to_wide()
# colnames(log_solomon345)
# 
# 
# KinEng_eq=log_solomon345[,2]%>%st_extract_spectra(alpha=0.05,del_num = 50) # stationary after 50 data points 
# print("This is Kinetic energy data")
# 
# PotEng_eq=log_solomon345[,3]%>%st_extract_spectra(alpha=0.05,del_num = 50) # stationary 
# print("This is Potential energy data")
# Press_eq=log_solomon345[,4]%>%st_extract_spectra(alpha=0.05,del_num = 50)# stationary 
# print("This is pressure data.")
# beep(8)
# Temp_eq=log_solomon345[,5]%>%st_extract_spectra(alpha=0.05,del_num = 50) # stationary after 50 data points 
# print("This is temperature data.")
# beep(8)
# TotEng_eq=log_solomon345[,6]%>%st_extract_spectra(alpha=0.05,del_num = 50) # stationary after 550 data points 
# print("This is total energy data")
# beep(8)

############################################################################################