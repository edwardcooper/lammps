## MSD g2 calculation for scaled positions a.k.a using xs,ys,zs,ix,iy,iz. 

### First define a function to do the calculation for one temeprature. 
MSD_scaled_g2_one_temp_Tmean=function(path="~/Dropbox/lammps/PS_20/atom300",filename="atom.300_1"
                                      ,molecule_atoms=645,num_mol=40
                                      ,atom_type=1:6,atom_type_mass=c(12.011,1.0079,12.011,12.011,12.011,1.0079),core_num=4){
  
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
  # calculate total atom number and timestep.
  
  tot_atom_num=atom.300.1_fread[,.(atom.id)]%>%max()
  paste("The number of atoms:",tot_atom_num,sep="")%>%message()
  
  timestep=dim(atom.300.1_fread)[1]/tot_atom_num
  
  paste("The number of timestep:",timestep,sep="")%>%message()
  
  ########################################################################
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
  # select only part of the data to release some memory. 
  atom.300.1_fread=atom.300.1_fread[,.(atom.id,type,time_step,xu,yu,zu)]
  
  
  ##########################################################################################################################################
  ##################################################################################
  # set atom.id and timestep both as the key for the data so that the data is sorted first by atom.id then by time_step
  atom.300.1_fread%>%setkey(atom.id,time_step)
  message("The data is sorted by atom.id and time_step")
  
  ##################################################################################
  
  
  #########################################################################################################
  # add mol variable to original data frame
  # Define a function to assign mol number according to the molecule size
  mol_id_gen=function(atom_id,molecule_atoms){
    # the minimum number for mol id is 0. 
    # minus 1 from the atom.id to get the last atom of each molecule to be on the correct molecule.
    return(floor((atom_id-1)/molecule_atoms))
  }
  
  timeRecordB()
  atom.300.1_fread[,mol:=mol_id_gen(atom_id = atom.id,molecule_atoms=molecule_atoms)]
  timeRecordB(output_message = "add mol variable data.table")
  #########################################################################################################
  # add mass variable to original data frame  
  # define a function to add mass to the dataset.
  
  decode = function(x, search, replace, default = NULL) {
    
    # build a nested ifelse function by recursion
    decode.fun <- function(search, replace, default = NULL)
      if (length(search) == 0) {
        function(x) if (is.null(default)) x else rep(default, length(x))
      } else {
        function(x) ifelse(x == search[1], replace[1],
                           decode.fun(tail(search, -1),
                                      tail(replace, -1),
                                      default)(x))
      }
    
    return(decode.fun(search, replace, default)(x))
  }
  
  timeRecordB()
  atom.300.1_fread[,mass:=decode(x = type, search = atom_type
                                 , replace = atom_type_mass)]
  timeRecordB(output_message = "add mass variable data.table")
  gc()
  
  #########################################################################################################
  
 
  
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
  
  #================================= This function assumes that the dataset contains a mol column. ==================================# 
  foreach_datas_lags_para_msd=function(data,lag_vec,mol_id,mol_min,core_num=core_num){
    registerDoParallel(cores = core_num)
    para_foreach_result=foreach(i=seq_along(mol_id),.combine = rbind)%dopar%{
      mol_center_of_mass=data[mol==mol_id[i],]%>%vapply_lags_msd(lag_vec = lag_vec)
      
      # if(i%%100==0 ){paste0("atom number",i,"/",length(mol_id))%>%message()} # echo the progress since I am a control freak. 
      return(mol_center_of_mass)
    }
    # getDoParWorkers()%>%print()
    stopImplicitCluster()
    gc()
    return(para_foreach_result)
  }
  foreach_datas_lags_para_msd=compiler::cmpfun(foreach_datas_lags_para_msd)
  
  # atom_id_1[,.(diff(xu,lag=2),diff(yu,lag=2),diff(zu,lag=2))]%>%head()%>%.^2%>%sum()
  
  ####################################################################################################################
  
  MSD_g2_matrix=function(data,timestep,num_mol,totmass,mol_min,core_num){
    # calculate the center of mass for every timestep and every mol 
    center_of_mass=data[,.(xu=sum(xu*mass)/totmass,yu=sum(yu*mass)/totmass,zu=sum(zu*mass)/totmass),by=.(mol,time_step)]
    # reorder the data first by mol then by the timestep
    center_of_mass%>%setkey(mol,time_step)
    # calculate the time average for MSD for each mol. 
    MSD_g2_empty_matrix=center_of_mass[,.(mol,xu,yu,zu)]%>%foreach_datas_lags_para_msd(lag_vec = 1:timestep,mol_id = mol_min:num_mol,core_num = core_num)
    
    return(MSD_g2_empty_matrix)
  } 
  
  #########################################################################################################
  # mass of molecule 1
  totmass= atom.300.1_fread[time_step==1,]%>%.[mol==1,.(mass)]%>%sum
  # calculate mol_min number
  mol_min=atom.300.1_fread[time_step==1,.(mol)]%>%min()

  # do it for every mol and combine the results by rbind
  
  timeRecordB()
  
  MSD.all.matrix=atom.300.1_fread%>%MSD_g2_matrix(timestep=timestep,num_mol=num_mol,totmass=totmass,core_num=core_num)
  
  gc()                                                                     
  
  timeRecordB(output_message = paste("cores:",core_num,"total mols:",num_mol, "MSD matrix calculation"))
  
  timeRecordR()%>%filter(output_message!="None")
  ########################
  
  # calculate non-gaussian parameter. 
  timeRecordB()
  NGP.600=(0.6)*colMeans(MSD.all.matrix^2)/(colMeans(MSD.all.matrix))^2-1
  timeRecordB(output_message = "NGP calculation")
  # write the NGP into a txt file.
  timeRecordB()
  NGP.600%>%write.table(file="NGP.g2.Tmean.1.txt", sep=",")
  timeRecordB(output_message = "NGP Write")
  
  # calculation the MSD averagd over all atoms and write it to a txt file.
  timeRecordB()
  MSD.all.matrix%>%colMeans%>%write.table(file="MSD.g2.colmean.Tmean.1.txt", sep=",")
  timeRecordB(output_message = "MSD average Write")
  
  
  gc()
  return(timeRecordR(ignore=0.1)%>%filter(output_message!="None")%>%select(output_message,run_time) )
  
}


# example use. 
# load the timeRecord functions from my github account.
source("https://raw.githubusercontent.com/edwardcooper/mlmodel_select/master/timeRecord_functions.R")
MSD_scaled_g2_one_temp_Tmean(path="~/Dropbox/lammps/PS_20/atom300",filename="atom.300_1"
                             ,molecule_atoms=645,num_mol=40,atom_type=1:6,atom_type_mass=c(12.011,1.0079,12.011,12.011,12.011,1.0079)
                            ,core_num=4)
#  stopping a parallel computation could only be done from htop. Not from the Rstudio.


## echo the current calculation and percentage of entire calculation.  
MSD_scaled_g2_Tmean=function(Path="~/Dropbox/lammps/",polymer="PS_20",temperatures=seq(200,520,by=20)
                             ,molecule_atoms=645,num_mol=40,atom_type=1:6,atom_type_mass=c(12.011,1.0079,12.011,12.011,12.011,1.0079)
                             ,core_num=4){
  # load the timeRecord functions from my github account.
  source("https://raw.githubusercontent.com/edwardcooper/mlmodel_select/master/timeRecord_functions.R")
  library(magrittr)
  # the loop to calculate the same thing in all temepratures defined above. 
  for (i in seq_along(temperatures)){
    # echo beginning of calculation
    paste("Begin calculation of temperature:",temperatures[i],sep="")%>%message
    
    # set correct path for the data file
    path=paste(Path,"/", polymer,"/atom",temperatures[i], sep='')
    # find the correct file to read and calculate.
    filename=paste("atom.",temperatures[i],"_1",sep="")
    
    # calculation for MSD
    MSD_scaled_g2_one_temp_Tmean(path=path,filename =filename
                                 ,molecule_atoms=molecule_atoms,num_mol=num_mol,
                                 atom_type=atom_type,atom_type_mass=atom_type_mass,polymer=polymer,core_num=core_num )
    
    # echo end of calculation
    paste("End calculation of temperature:",temperatures[i],sep="")%>%message
    paste(i,"/",length(temperatures))%>%message
    gc()
  }
  
  return( timeRecordR(ignore=0.1)%>%filter(output_message!="None")%>%select(output_message,run_time) )
}

 # example use
 # MSD_scaled_g2_Tmean(Path="~/Dropbox/lammps/",polymer="PS_20",temperatures=seq(200,520,by=20)
 #                    ,molecule_atoms=645,num_mol=40,atom_type=1:6,atom_type_mass=c(12.011,1.0079,12.011,12.011,12.011,1.0079)
 #                    ,core_num=4)

