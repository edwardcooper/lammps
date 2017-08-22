# Stationarity test on total energy, kinetic energy, potential energy, pressure and temperature data.
# This is for the purpose to see if the system is already in equilibrium. 


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
  dataset=cbind(timestep=1:dim(dataset)[1],dataset)
  gc()
  return(dataset)
  
}
log_data=read.table("log_npt.csv")
log_data=log_data%>%long_to_wide()

colnames(log)

# After loading packages and data, we will need to do hypthesis tests here.
###########################################################################



#var_num=dim(log_data)[2]-1

test=fractal::stationarity(log_data[,2],n.block = 10)
summary(test)

tseries::adf.test(log_data[,2])
tseries::kpss.test(log_data[,2])
tseries::pp.test(log_data[,2])

stats::PP.test(log_data[,2])

aTSA::adf.test(log_data[,2])
aTSA::pp.test(log_data[,2])
aTSA::kpss.test(log_data[,2])


# This takes a very very long time. 
# urca::ur.za(log[,2])

# It would be best avoid this method since it has some length restrictions.
# library(locits)
# # The length of time series must be a number that is 2^n. 
# # Use log2(length) to find a working length for this algorithm.
# test2=hwtos2(Kineng[1:2^15])
# 
# test2

