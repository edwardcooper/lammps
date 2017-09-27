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