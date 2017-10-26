## Write a function to calculate MSD averaged over all atoms for all temperatures. 

### First define a function to do the calculation for one temeprature. 
MSD_g0_one_temp=function(path="~/Dropbox/lammps/PMMA_big/atom600",filename="atom.600_1.txt",timestep=5001){
  # load the timeRecord functions from my github account.
  source("https://raw.githubusercontent.com/edwardcooper/mlmodel_select/master/timeRecord_functions.R")
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
  atom.600.1_fread=fread(input=filename, sep=" ", stringsAsFactors = FALSE, fill=TRUE
                         #,col.names = c("atom-id","type","mol","xu","yu","zu")
                         ,colClasses = c("numeric","numeric","numeric","numeric","numeric","numeric","character","character")
                         ,na.strings = c("ITEM:","TIMESTEP","NUMBER","OF","ATOMS","BOX","BOUNDS", "pp","id","type","mol","xu","yu","zu")
                         
  )
  timeRecordB(output_message = "Load in data fread")
  # select the non-NA columns
  atom.600.1_fread=atom.600.1_fread[,.(V1,V2,V3,V4,V5,V6)]
  # clear the rows that have NA values, V6 contains the most NAs. 
  atom.600.1_fread=atom.600.1_fread[complete.cases(atom.600.1_fread[,V6])]
  # set the column names. 
  colnames(atom.600.1_fread)=c("atom.id","type","mol","xu","yu","zu")
  # define the total number of atoms. 
  tot_atom_num=atom.600.1_fread[,"atom.id"]%>%max()
  # Define MSD function
  MSD=function(data){
    (data$xu-data$xu[1])^2+ (data$yu-data$yu[1])^2+ (data$zu-data$zu[1])^2
  }
  
  # Define a function to combine the results from MSD calculation
  MSD_matrix=function(data,timestep,tot_atoms){
    MSD_empty_matrix=matrix(NA,nrow=tot_atoms,ncol=timestep)
    for(i in 1:tot_atoms){
      MSD_empty_matrix[i,]=data%>%.[.(i),.(xu,yu,zu)]%>%MSD
    }
    return(MSD_empty_matrix)
  }
  # record the time. 
  timeRecordB()
  # set the key for faster search and shorter codes. 
  setkey(atom.600.1_fread,atom.id)
  # Do the MSD calculation for every atom. 
  MSD.all.matrix.600=atom.600.1_fread%>%MSD_matrix(timestep = timestep,tot_atoms = tot_atom_num)
  # record the time again.
  timeRecordB(output_message = "MSD matrix calculation for with replacement")
  
  # calculate non-gaussian parameter. 
  timeRecordB()
  NGP.600=(0.6)*colMeans(MSD.all.matrix.600^2)/(colMeans(MSD.all.matrix.600))^2-1
  timeRecordB(output_message = "NGP calculation")
  # write the NGP into a txt file
  timeRecordB()
  NGP.600%>%write.table(file="NGP.g0.1.txt", sep=",")
  timeRecordB(output_message = "NGP Write")
  
  # calculation the MSD averagd over all atoms and write it to a txt file.
  timeRecordB()
  MSD.all.matrix.600%>%colMeans%>%write.table(file="MSD.g0.colmean.1.txt", sep=",")
  timeRecordB(output_message = "MSD average Write")
  
  return(timeRecordR(ignore=0.1)%>%filter(output_message!="None")%>%select(output_message,run_time)
)
}






# test usage()
# MSD_g0_one_temp(path="~/Dropbox/lammps/PMMA_big/atom300",filename="atom.300_1.txt",timestep=5001)




## echo the current calculation and percentage of entire calculation.  
MSD_g0=function(Path="~/Dropbox/lammps/",polymer="PMMA_big",temperatures=seq(300,600,by=20),timestep=5001){
  library(magrittr)
  # the loop to calculate the same thing in all temepratures defined above. 
  for (i in seq_along(temperatures)){
    # echo beginning of calculation
    paste("Begin calculation of temperature:",temperatures[i],sep="")%>%message
    
    # set correct path for the data file
    path=paste(Path,"/", polymer,"/atom",temperatures[i], sep='')
    # find the correct file to read and calculate.
    filename=paste("atom.",temperatures[i],"_1.txt",sep="")
    
    # calculation for MSD
    MSD_g0_one_temp(path=path,filename =filename ,timestep=timestep)
    
    # echo end of calculation
    paste("End calculation of temperature:",temperatures[i],sep="")%>%message
    paste(i,"/",length(temperatures))%>%message
  }
   
  return( timeRecordR(ignore=0.1)%>%filter(output_message!="None")%>%select(output_message,run_time) )
}

# example use 
# MSD_g0()