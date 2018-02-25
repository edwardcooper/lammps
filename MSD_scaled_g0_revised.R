## MSD g0 calculation for scaled positions a.k.a using xs,ys,zs,ix,iy,iz. 


### First define a function to do the calculation for one temeprature. 
MSD_scaled_g0_one_temp=function(path="~/Dropbox/lammps/PMMA_long/atom300",filename="atom.300_long2"){
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
  atom.300.1_fread=fread(input=filename, sep=" ", stringsAsFactors = FALSE, fill=TRUE
                         #,col.names = c("atom-id","type","xs","ys","zs","ix","iy","iz")
                         ,colClasses = c("numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","character","character")
                         ,na.strings = c("ITEM:","TIMESTEP","NUMBER","OF","ATOMS","BOX","BOUNDS", "pp","id","type","mol","xs","ys","zs","ix","iy","iz")
  )
  ## first give the bounds of box on the x,y,z direction. 
  xhigh=atom.300.1_fread[6,V2]
  xlow=atom.300.1_fread[6,V1]
  yhigh=atom.300.1_fread[7,V2]
  ylow=atom.300.1_fread[7,V1]
  zhigh=atom.300.1_fread[8,V2]
  zlow=atom.300.1_fread[8,V1]
  
  paste("The box bounds are:\n","x:",xlow,xhigh,"\n","y:",ylow,yhigh,"\n","z:",zlow,zhigh,"\n")%>%message
  
  
  timeRecordB(output_message = "Load in data fread")
  # select the non-NA columns
  atom.300.1_fread=atom.300.1_fread[,.(V1,V2,V3,V4,V5,V6,V7,V8)]
  # clear the rows that have NA values, V6 contains the most NAs. 
  atom.300.1_fread=atom.300.1_fread[complete.cases(atom.300.1_fread[,V6])]
  colnames(atom.300.1_fread)=c("atom.id","type","xs","ys","zs","ix","iy","iz")
  
  
  
  ########################################################################
  
  # Define MSD function
  MSD=function(data){
    (data$xu-data$xu[1])^2+ (data$yu-data$yu[1])^2+ (data$zu-data$zu[1])^2
  }
  
  
  tot_atom_num=atom.300.1_fread[,.(atom.id)]%>%max()
  
  timestep=dim(atom.300.1_fread)[1]/tot_atom_num
  
  paste("The number of timestep:",timestep,sep="")%>%message()
  ##################################################################################
  ### develop along this line the average MSD between all atoms.####################
  ##################################################################################
  
  # add time variable for easier data handling and safe guarding any unwanted ordering to change the ordering of time in the data. 
  timeRecordB()
  atom.300.1_fread[,time_step:=seq(1,timestep,by=1)%>%rep(times=tot_atom_num)%>%sort()]
  timeRecordB(output_message = "add time variable data.table")
  
  
  # calculate the xu,yu,zu unscaled position 
  
  ## first give the bounds of box on the x,y,z direction. 
  
  
  
  xlength=xhigh-xlow
  ylength=yhigh-ylow
  zlength=zhigh-zlow
  
  
  atom.300.1_fread[,`:=`(xu=xs*xlength+xlow+ix*xlength,yu=ys*ylength+ylow+iy*ylength,zu=zs*zlength+zlow+iz*zlength)]
  
  #####################################################################
  MSD_matrix=function(data,timestep,tot_atoms){
    MSD_empty_matrix=matrix(NA,nrow=tot_atoms,ncol=timestep)
    for(i in 1:tot_atoms){
      MSD_empty_matrix[i,]=data[atom.id==i,.(xu,yu,zu)]%>%MSD
    }
    return(MSD_empty_matrix)
  }
  timeRecordB()
  
  
  MSD.all.matrix.300=atom.300.1_fread%>%MSD_matrix(timestep = timestep,tot_atoms = tot_atom_num)
  timeRecordB(output_message = "MSD matrix calculation for with replacement")
  
  
  
  timeRecordB()
  NGP.300=(0.6)*colMeans(MSD.all.matrix.300^2)/(colMeans(MSD.all.matrix.300))^2-1
  timeRecordB(output_message = "NGP calculation")
  
  timeRecordB()
  NGP.300%>%write.table(file="NGP.g0.1.txt", sep=",")
  timeRecordB(output_message = "NGP Write")
  
  
  timeRecordB()
  MSD.all.matrix.300%>%colMeans%>%write.table(file="MSD.g0.colmean.1.txt", sep=",")
  timeRecordB(output_message = "MSD average over all atom calculation and write")
  
  
  return(timeRecordR(ignore=0.1)%>%filter(output_message!="None")%>%select(output_message,run_time)
  )
}


# example use 
# MSD_scaled_g0_one_temp(path="~/Dropbox/lammps/PS_20/atom300",filename="atom.300_1"
#                       )


## echo the current calculation and percentage of entire calculation.  
MSD_scaled_g0=function(Path="~/Dropbox/lammps/",polymer="PMMA_long",temperatures=seq(300,620,by=20)){
  library(magrittr)
  # the loop to calculate the same thing in all temepratures defined above. 
  for (i in seq_along(temperatures)){
    # echo beginning of calculation
    paste("Begin calculation of temperature:",temperatures[i],sep="")%>%message
    
    # set correct path for the data file
    path=paste(Path,"/", polymer,"/atom",temperatures[i], sep='')
    # find the correct file to read and calculate.
    filename=paste("atom.",temperatures[i],"_long2",sep="")
    
    # calculation for MSD
    MSD_scaled_g0_one_temp(path=path,filename =filename,xhigh=xhigh,xlow=xlow,yhigh=yhigh,ylow=ylow,zhigh=zhigh,zlow=zlow )
    
    # echo end of calculation
    paste("End calculation of temperature:",temperatures[i],sep="")%>%message
    paste(i,"/",length(temperatures))%>%message
    gc()
  }
  
  return( timeRecordR(ignore=0.1)%>%filter(output_message!="None")%>%select(output_message,run_time) )
}

# example use
# MSD_scaled_g0(Path="~/Dropbox/lammps/",polymer="PMMA_long",temperatures=seq(300,620,by=20) )

