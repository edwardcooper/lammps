
### First define a function to do the calculation for one temeprature. 
MSD_unscaled_g0_one_temp_Tmean=function(path="~/Dropbox/lammps/PMMA_big/atom300",filename="atom.300_1.txt",core_num=10){
  # load the timeRecord functions from my github account.
  setwd(path)
  # record the time
  timeRecordB()
  # load all library
  library(data.table)
  library(foreach)
  library(dplyr)
  library(matrixStats)
  library(magrittr)
  ### Load the data with fread
  atom.300.1_fread=fread(input=filename, sep=" ", stringsAsFactors = FALSE, fill=TRUE
                         #,col.names = c("atom-id","type","mol","xu","yu","zu")
                         ,colClasses = c("numeric","numeric","numeric","numeric","numeric","numeric","character","character")
                         ,na.strings = c("ITEM:","TIMESTEP","NUMBER","OF","ATOMS","BOX","BOUNDS", "pp","id","type","mol","xu","yu","zu")
                         
  )
  timeRecordB(output_message = "Load in data fread")
  # select the non-NA columns
  atom.300.1_fread=atom.300.1_fread[,.(V1,V2,V3,V4,V5,V6)]
  # clear the rows that have NA values, V6 contains the most NAs. 
  atom.300.1_fread=atom.300.1_fread[complete.cases(atom.300.1_fread[,V6])]
  # set the column names. 
  colnames(atom.300.1_fread)=c("atom.id","type","mol","xu","yu","zu")
  
  
  timeRecordB(output_message = "Load in data fread")
  
  
  
  ########################################################################
  # calculate total atom number and timestep.
  
  tot_atom_num=atom.300.1_fread[,.(atom.id)]%>%max()
  
  timestep=dim(atom.300.1_fread)[1]/tot_atom_num
  
  paste("The number of timestep:",timestep,sep="")%>%message()
  
  ########################################################################
  # add time variable for easier data handling and safe guarding any unwanted ordering to change the ordering of time in the data. 
  timeRecordB()
  atom.300.1_fread[,time_step:=seq(1,timestep,by=1)%>%rep(times=tot_atom_num)%>%sort()]
  timeRecordB(output_message = "add time variable data.table")
  
  ##########################################################################################################################################
  
  #########################################################################################################
  # define functions to calculate time average 
  
  library(foreach)
  library(doParallel)
  library(matrixStats)
  library(magrittr)
  
  # calculate the lag difference, then squared all the values, calculate the column mean for sqaured difference, then sum the difference in x,y,z directions.
  coldiff_msd=function(data,lag_val){
    return(data%>%as.matrix%>%colDiffs(lag=lag_val)%>%.^2%>%colMeans%>%sum)
  }
  coldiff_msd=compiler::cmpfun(coldiff_msd)
  
  # calculate the msd for a vector of lag values
  vapply_lags_msd=function(data,lag_vec){
    vapply(lag_vec,function(lag_val){coldiff_msd(data,lag_val = lag_val)},FUN.VALUE = numeric(1))
  }
  vapply_lags_msd=compiler::cmpfun(vapply_lags_msd)
  
  # calculate lagged msd for different atom.id. 
  # This the most outerloop per se. So that we could use foreach to do things in parallel. 
  # This function assumes that the dataset contains a atom.id column. 
  foreach_datas_lags_para_msd=function(data,lag_vec,atom_id,core_count=2){
    registerDoParallel(cores = core_count)
    para_foreach_result=foreach(i=seq_along(atom_id),.combine = rbind)%dopar%{
      data[atom.id==atom_id[i],]%>%vapply_lags_msd(lag_vec = lag_vec)
      # if(i%%1000==0 ){paste0("atom number",i,"/",tot_atom_num)} # if you want to add this, do not forget a return value for the calculation or nothing will return.
    }
    # getDoParWorkers()%>%print()
    stopImplicitCluster()
    gc()
    return(para_foreach_result)
  }
  foreach_datas_lags_para_msd=compiler::cmpfun(foreach_datas_lags_para_msd)
  # all the functions defined and compiled for calculating the time average. 
  ################################################################################################
  
  
  ###################################################################################################
  # start the main computation part!
  
  timeRecordB()
  # do it for every atom.id and combine the results by rbind
  MSD.all.matrix=atom.300.1_fread[,.(atom.id,xu,yu,zu)]%>%foreach_datas_lags_para_msd(lag_vec = 1:timestep,atom_id = 1:tot_atom_num,core_count = core_num)
  
  timeRecordB(output_message = paste("cores:",core_num,"total atoms:",tot_atom_num, "MSD matrix calculation"))
  
  
  
  # main computation part ends. 
  ###################################################################################################
  
  
  
  
  # calculate non-gaussian parameter. 
  timeRecordB()
  NGP.300=(0.6)*colMeans(MSD.all.matrix^2)/(colMeans(MSD.all.matrix))^2-1
  timeRecordB(output_message = "NGP calculation")
  # write the NGP into a txt file
  timeRecordB()
  NGP.300%>%write.table(file="NGP.g0.Tmean.1.txt", sep=",")
  timeRecordB(output_message = "NGP Write")
  
  # calculation the MSD averagd over all atoms and write it to a txt file.
  timeRecordB()
  MSD.all.matrix%>%colMeans%>%write.table(file="MSD.g0.colmean.Tmean.1.txt", sep=",")
  timeRecordB(output_message = "MSD average Write")
  
  
  return(timeRecordR(ignore=0.1)%>%filter(output_message!="None")%>%select(output_message,run_time) )
  
}


# example use 
# source("https://raw.githubusercontent.com/edwardcooper/mlmodel_select/master/timeRecord_functions.R")
# MSD_unscaled_g0_one_temp_Tmean(path="~/Dropbox/lammps/PMMA_big/atom300",filename="atom.300_1.txt",core_num = 10)
# stopping a parallel computation could only be done from htop. Not from the Rstudio.


## echo the current calculation and percentage of entire calculation.  
MSD_unscaled_g0_Tmean=function(Path="~/Dropbox/lammps/",polymer="PS_20",temperatures=seq(200,520,by=20),core_num=10){
  library(magrittr)
  source("https://raw.githubusercontent.com/edwardcooper/mlmodel_select/master/timeRecord_functions.R")
  # the loop to calculate the same thing in all temepratures defined above. 
  for (i in seq_along(temperatures)){
    # echo beginning of calculation
    paste("Begin calculation of temperature:",temperatures[i],sep="")%>%message
    
    # set correct path for the data file
    path=paste(Path,"/", polymer,"/atom",temperatures[i], sep='')
    # find the correct file to read and calculate.
    filename=paste("atom.",temperatures[i],"_1.txt",sep="")
    
    # calculation for MSD
    MSD_unscaled_g0_one_temp_Tmean(path=path,filename =filename,core_num = core_num )
    
    # echo end of calculation
    paste("End calculation of temperature:",temperatures[i],sep="")%>%message
    paste(i,"/",length(temperatures))%>%message
    gc()
  }
  
  return( timeRecordR(ignore=0.1)%>%filter(output_message!="None")%>%select(output_message,run_time) )
}

# example use
MSD_unscaled_g0_Tmean(Path="~/Dropbox/lammps/",polymer="PMMA_big",temperatures=seq(300,600,by=20),core_num = 10)
