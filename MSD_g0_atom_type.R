## Write a function to calculate MSD averaged over all atoms for a certain atom type.


## First define a function to do the calculation for one temeprature. 
MSD_g0_atom_type_one_temp=function(path="~/Dropbox/lammps/PMMA_big/atom600",filename="atom.600_1.txt"){
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
  #################################################################################
  
  # calculate the timestep from the data ( divide the total number of rows by the number of atoms)
  timestep=dim(atom.600.1_fread)[1]/max(atom.600.1_fread[,atom.id])
  # calculate total number of atoms assuming that smallest atom.id is 1. 
  tot_atom_num=max(atom.600.1_fread[,atom.id])
  #################################################################################
  
  # add time variable for easier data handling and safe guarding any unwanted ordering to change the ordering of time in the data. 
  timeRecordB()
  atom.600.1_fread[,time_step:=seq(1,timestep,by=1)%>%rep(times=tot_atom_num)%>%sort()]
  timeRecordB(output_message = "add time variable data.table")
  
  #################################################################################
  # Define MSD function
  MSD=function(data){
    (data$xu-data$xu[1])^2+ (data$yu-data$yu[1])^2+ (data$zu-data$zu[1])^2
  }
  
  # Define a function for calculating MSD, NGP for all atom types
  atom_type_g0_MSD_NGP=function(data){
    # calculate parameters based on the data
    typemin=data[,type]%>%max
    typemax=data[,type]%>%min
    #number of atoms
    tot_atom_num=max(data[,atom.id])
    #timestep
    timestep=timestep=dim(atom.600.1_fread)[1]/max(atom.600.1_fread[,atom.id])
    
    # MSD calculation for every atom type
    #######################################################################################################
    
    for( n in typemin:typemax){
      
      # find the atom.id for type 1 atoms.
      atom_id_for_different_type=data[time_step==1,]%>%.[type==n,atom.id]
      
      #####################################################################
      MSD_matrix=function(data,timestep,tot_atoms){
        MSD_empty_matrix=matrix(NA,nrow=tot_atoms,ncol=timestep)
        for(i in 1:tot_atoms){
          MSD_empty_matrix[i,]=data%>%.[atom.id==atom_id_for_different_type[i],.(xu,yu,zu)]%>%MSD
        }
        return(MSD_empty_matrix)
      }
      timeRecordB()
      
      
      atom.type.MSD.matrix=atom.600.1_fread%>%MSD_matrix(timestep = timestep,tot_atoms = length(atom_id_for_different_type))
      out_message=paste("MSD matrix calculation for with replacement for atoms of type ",n,sep="")
      timeRecordB(output_message =out_message )
      paste("The total atoms for type ",n," is ",length(atom_id_for_different_type))%>%message()
      # calculate colmean, NGP and write it onto disk.
      #####################################################################
      MSD.colmean.name=paste("atom.type.",n,".g0_MSD.colmean.txt", sep='')
      atom.type.MSD.matrix%>% colMeans%>%write.table(file=MSD.colmean.name,sep=",")
      
      NGP.name=paste("atom.type.",n,".g0_NGP.txt", sep='')
      NGP=(0.6)*colMeans(atom.type.MSD.matrix^2)/(colMeans(atom.type.MSD.matrix))^2-1
      NGP%>%write.table(file=NGP.name,sep=",")
    }
    
  }
  

  # do the calculation for all atom types
  atom_type_g0_MSD_NGP(data=atom.600.1_fread)

  
  return(timeRecordR(ignore=0.1)%>%filter(output_message!="None")%>%select(output_message,run_time)
  )
}

# example use 
# MSD_g0_atom_type_one_temp(path="~/Dropbox/lammps/PMMA_big/atom500",filename="atom.500_1.txt")


## echo the current calculation and percentage of entire calculation.  
MSD_g0_atom_type=function(Path="~/Dropbox/lammps/",polymer="PMMA_big",temperatures=seq(300,600,by=20)){
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
    MSD_g0_atom_type_one_temp(path=path,filename =filename)
    
    # echo end of calculation
    paste("End calculation of temperature:",temperatures[i],sep="")%>%message
    paste(i,"/",length(temperatures))%>%message
    gc()
  }
  
  return( timeRecordR(ignore=0.1)%>%filter(output_message!="None")%>%select(output_message,run_time) )
}

# example use 
MSD_g0_atom_type(Path="~/Dropbox/lammps/",polymer="PMMA_big",temperatures=seq(300,600,by=20))
