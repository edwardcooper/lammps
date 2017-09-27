########################################################################################
# The spread function from tidyr gives me an error "Error: C stack usage  25955701 is too close to the limit".
# I have not been able to debug the problem. So I will have to write my own R function to do it.



# This function mimics the usage of spread function,
# which convert a dataset with a long format into a wide format.

# This function only works with a two-column-data case
# where the first column is the name of the variable, second column is the value of the variable.

# for more details concerning what is spread and gather function in tidyr, please refer to http://r4ds.had.co.nz/tidy-data.html#spreading-and-gathering


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

# log=read.table("log.csv")
# # This is an example of how to use this function
# log=log%>%long_to_wide()
# 
# head(log)
# colnames(log)

# I will try to extend its functionality if necessary in the future.

#########################################################################
# read in dataset from the log file to examine the results 



log_plot=function(data="NULL",file="NULL"){
  if(file!="NULL"){
    # read in the data
    log=read.table(file)
    print("Finished loading data.")
    log=log%>%long_to_wide()
    print("Finished converting data to wide format.")
  }else if(data!="NULL"){
    log=data
  }else{
    print("Enter either a filename or a dataset")
  }
  
  
  library(ggplot2)
  
  p1=log%>%ggplot()+geom_line(aes(x=timestep,y=Temp),color="blue")
  p2=log%>%ggplot()+geom_line(aes(x=timestep,y=Press),color="green")
  
  p3=log%>%ggplot()+geom_line(aes(x=timestep,y=KinEng),color="orange")
  p4=log%>%ggplot()+geom_line(aes(x=timestep,y=PotEng),color="red")
  p5=log%>%ggplot()+geom_line(aes(x=timestep,y=TotEng),color="black")
  
  print("Finished plotting separate graphs")
  library(grid)
  grid.newpage()
  plot1=grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "first"))
  
  
  grid.newpage()
  plot2=grid.draw(rbind(ggplotGrob(p3), ggplotGrob(p4),ggplotGrob(p5), size = "first"))
  
  gc()
  return(print("Finished plotting"))
}




# log_plot("log_mini.csv")
# log_plot("log_npt.csv")
# log_plot("log_npt_noS.csv")
# log_plot("log320_solomon.csv")
# log_plot("log345_solomon.csv")
