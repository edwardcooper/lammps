## MSD g1 calculation for unscaled positions a.k.a using xs,ys,zs,ix,iy,iz. 

### First define a function to do the calculation for one temeprature. 
MSD_unscaled_g1_one_temp_Tmean=function(path="~/Dropbox/lammps/PMMA_big/atom300"
                                        ,filename="atom.300_1.txt"              
                                        ,polymer="PMMA_big"
                                        ,num_mol=64
                                        ,molecule_atoms=602
                                        ,molecule_monomers=40
                                        ,monomer_atoms=15
                                        ,atom_type=1:10
                                        ,atom_type_mass=c(1.0079,12.011,12.011,12.011,15.9999,15.9999,12.011,12.011,1.0079,12.011),core_num=4){
  
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
  paste("The number of atoms:",tot_atom_num,sep="")%>%message()
  
  timestep=dim(atom.300.1_fread)[1]/tot_atom_num
  
  paste("The number of timestep:",timestep,sep="")%>%message()
  
  ########################################################################
  # add time variable for easier data handling and safe guarding any unwanted ordering to change the ordering of time in the data. 
  timeRecordB()
  atom.300.1_fread[,time_step:=seq(1,timestep,by=1)%>%rep(times=tot_atom_num)%>%sort()]
  timeRecordB(output_message = "add time variable data.table")
  
  

  
  
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
  
  ##########
  # define a function to generate correct monomer id depending on the polymer option in the function.
  if(polymer %in% c("PMMA_big","PMMA_long")){
    
    
    monomer_gen=function(atom_id,molecule_atoms,molecule_monomers,monomer_atoms,edge_atoms=c(0,1)){
      # if it is on the either end of polymer, the monomer is defined as monomer 0. 
      # if the atom is not on the end, then first calculate the molecule number and multiply it by 40 since there is 40 monomers in each molecule.
      # then add the the monomer number it has in this molecule. Need to deduct the first atom off the molecule, divided by the number of atoms in a monomer and get a integer
      # then you have the monomer id. 
      monomer.id=ifelse(atom_id%%molecule_atoms %in% edge_atoms
                        ,0
                        ,floor(atom_id/molecule_atoms)*molecule_monomers+ceiling((atom_id%%molecule_atoms-1)/monomer_atoms) )
      return(monomer.id)
    }
    
    
  }else if(polymer %in% c("PS","PS_20") ){
    
    monomer_gen=function(atom_id,molecule_atoms,molecule_monomers,monomer_atoms,edge_atoms=c(17,642,643,644,0)){
      
      monomer.id=ifelse((atom_id%%molecule_atoms )%in% edge_atoms
                        ,0
                        ,floor(atom_id/molecule_atoms)*molecule_monomers
                        +ceiling((atom_id%%molecule_atoms+ifelse(atom_id%%molecule_atoms>16,-1,0))/monomer_atoms)
      )
      return(monomer.id)
    }
    
  }
  
  ## add here more else if when there is more polymers. 
  ##########
  # add monomer id 
  timeRecordB()
  atom.300.1_fread[,monomer.id:=monomer_gen(atom_id=atom.id,molecule_atoms=molecule_atoms,molecule_monomers=molecule_monomers,monomer_atoms=monomer_atoms)]
  timeRecordB(output_message = "add monomer id variable data.table")
  gc()
  #########################################################################################################
  # mass of monomer 
  monomer_mass=atom.300.1_fread[time_step==1,]%>%.[monomer.id==1,.(mass)]%>%sum
  
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
  
  #================================= This function assumes that the dataset contains a monomer.id column. ==================================# 
  foreach_datas_lags_para_msd=function(data,lag_vec,monomer_id,core_num=core_num){
    registerDoParallel(cores = core_num)
    para_foreach_result=foreach(i=seq_along(monomer_id),.combine = rbind)%dopar%{
      monomer_center_of_mass=data[monomer.id==monomer_id[i],]%>%vapply_lags_msd(lag_vec = lag_vec)
      
      # if(i%%100==0 ){paste0("atom number",i,"/",length(monomer_id))%>%message()} # echo the progress since I am a control freak. 
      return(monomer_center_of_mass)
    }
    # getDoParWorkers()%>%print()
    stopImplicitCluster()
    gc()
    return(para_foreach_result)
  }
  foreach_datas_lags_para_msd=compiler::cmpfun(foreach_datas_lags_para_msd)
  
  # atom_id_1[,.(diff(xu,lag=2),diff(yu,lag=2),diff(zu,lag=2))]%>%head()%>%.^2%>%sum()
  
  ####################################################################################################################
  
  MSD_g1_matrix=function(data,timestep,num_monomer,monomer_mass,tot_monomer_num,core_num){
    # calculate the center of mass for every timestep and every monomer 
    center_of_mass=data[,.(xu=sum(xu*mass)/monomer_mass,yu=sum(yu*mass)/monomer_mass,zu=sum(zu*mass)/monomer_mass),by=.(monomer.id,time_step)]
    # reorder the data first by monomer.id then by the timestep
    center_of_mass%>%setkey(monomer.id,time_step)
    # calculate the time average for MSD for each monomer. 
    MSD_g1_empty_matrix=center_of_mass[,.(monomer.id,xu,yu,zu)]%>%foreach_datas_lags_para_msd(lag_vec = 1:timestep,monomer_id = 1:tot_monomer_num,core_num = core_num)
    
    return(MSD_g1_empty_matrix)
  } 
  
  #########################################################################################################
  # calculate the total monomer number. 
  tot_monomer_num=num_mol*molecule_monomers
  
  # do it for every monomer.id and combine the results by rbind
  
  timeRecordB()
  
  MSD.all.matrix=atom.300.1_fread%>%MSD_g1_matrix(timestep=timestep,num_monomer=num_mol*molecule_monomers,monomer_mass=monomer_mass
                                                  ,tot_monomer_num=tot_monomer_num,core_num=core_num)
  
  gc()                                                                     
  
  timeRecordB(output_message = paste("cores:",core_num,"total monomers:",tot_monomer_num, "MSD matrix calculation"))
  
  timeRecordR()%>%filter(output_message!="None")
  ########################
  
  # calculate non-gaussian parameter. 
  timeRecordB()
  NGP.600=(0.6)*colMeans(MSD.all.matrix^2)/(colMeans(MSD.all.matrix))^2-1
  timeRecordB(output_message = "NGP calculation")
  # write the NGP into a txt file.
  timeRecordB()
  NGP.600%>%write.table(file="NGP.g1.Tmean.1.txt", sep=",")
  timeRecordB(output_message = "NGP Write")
  
  # calculation the MSD averagd over all atoms and write it to a txt file.
  timeRecordB()
  MSD.all.matrix%>%colMeans%>%write.table(file="MSD.g1.colmean.Tmean.1.txt", sep=",")
  timeRecordB(output_message = "MSD average Write")
  
  
  
  return(timeRecordR(ignore=0.1)%>%filter(output_message!="None")%>%select(output_message,run_time) )
  
}


# example use. 
# load the timeRecord functions from my github account.
# source("https://raw.githubusercontent.com/edwardcooper/mlmodel_select/master/timeRecord_functions.R")
# MSD_unscaled_g1_one_temp_Tmean(path="~/Dropbox/lammps/PMMA_big/atom300"
#                                ,filename="atom.300_1.txt"              
#                                ,polymer="PMMA_big"
#                                ,num_mol=64
#                                ,molecule_atoms=602
#                                ,molecule_monomers=40
#                                ,monomer_atoms=15
#                                ,atom_type=1:10
#                                ,atom_type_mass=c(1.0079,12.011,12.011,12.011,15.9999,15.9999,12.011,12.011,1.0079,12.011),core_num=4)
#  stopping a parallel computation could only be done from htop. Not from the Rstudio.


## echo the current calculation and percentage of entire calculation.  
MSD_unscaled_g1_Tmean=function(Path="~/Dropbox/lammps/",polymer="PMMA_big",temperatures=seq(300,600,by=20)
                               ,num_mol=64
                               ,molecule_atoms=602
                               ,molecule_monomers=40
                               ,monomer_atoms=15
                               ,atom_type=1:10
                               ,atom_type_mass=c(1.0079,12.011,12.011,12.011,15.9999,15.9999,12.011,12.011,1.0079,12.011)
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
    filename=paste("atom.",temperatures[i],"_1.txt",sep="")
    
    # calculation for MSD
    MSD_unscaled_g1_one_temp_Tmean(path=path,filename =filename
                                 ,molecule_atoms=molecule_atoms,num_mol=num_mol,molecule_monomers=molecule_monomers
                                 ,monomer_atoms=monomer_atoms,atom_type=atom_type,atom_type_mass=atom_type_mass,polymer=polymer,core_num=core_num )
    
    # echo end of calculation
    paste("End calculation of temperature:",temperatures[i],sep="")%>%message
    paste(i,"/",length(temperatures))%>%message
    gc()
  }
  
  return( timeRecordR(ignore=0.1)%>%filter(output_message!="None")%>%select(output_message,run_time) )
}

# example use
# MSD_unscaled_g1_Tmean(Path="~/Dropbox/lammps/",polymer="PMMA_big",temperatures=seq(300,600,by=20)
#                       ,num_mol=64
#                       ,molecule_atoms=602
#                       ,molecule_monomers=40
#                       ,monomer_atoms=15
#                       ,atom_type=1:10
#                       ,atom_type_mass=c(1.0079,12.011,12.011,12.011,15.9999,15.9999,12.011,12.011,1.0079,12.011)
#                       ,core_num=4)
# 
# 
