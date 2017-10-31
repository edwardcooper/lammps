# MSD calculation for scaled position coordinates xs, ix .

### First define a function to do the calculation for one temeprature. 
MSD_scaled_g1_one_temp=function(path="~/Dropbox/lammps/PMMA_long/atom300",filename="atom.300_long2",polymer="PMMA_long"
                                ,num_mol=64,molecule_atoms=602,molecule_monomers=40 ,monomer_atoms=15
                                ,atom_type=1:10,atom_type_mass=c(1.0079,12.011,12.011,12.011,15.9999,15.9999,12.011,12.011,1.0079,12.011)
                                ,xhigh=72.3695
                                ,xlow=2.86271
                                ,yhigh=75.7675
                                ,ylow=-0.535315
                                ,zhigh=73.6408
                                ,zlow=1.59139
                                ){
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
  timeRecordB(output_message = "Load in data fread")
  # select the non-NA columns
  atom.300.1_fread=atom.300.1_fread[,.(V1,V2,V3,V4,V5,V6,V7,V8)]
  # clear the rows that have NA values, V6 contains the most NAs. 
  atom.300.1_fread=atom.300.1_fread[complete.cases(atom.300.1_fread[,V6])]
  #set the column names
  colnames(atom.300.1_fread)=c("atom.id","type","xs","ys","zs","ix","iy","iz")
  
  
  
  ########################################################################
  
  # Define MSD function
  MSD=function(data){
    (data$xu-data$xu[1])^2+ (data$yu-data$yu[1])^2+ (data$zu-data$zu[1])^2
  }
  #define a function for replacement 
  
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
  
  # Define a function to assign mol number according to the molecule size
  mol_id_gen=function(atom_id,molecule_atoms){
    # the minimum number for mol id is 0. 
    # minus 1 from the atom.id to get the last atom of each molecule to be on the correct molecule.
    return(floor((atom_id-1)/molecule_atoms))
  }
  
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
    
    
  }else if(polymer=="PS"){
    
    monomer_gen=function(atom_id,molecule_atoms,molecule_monomers,monomer_atoms,edge_atoms=c(17,642,643,644,0)){
      
      monomer.id=ifelse((atom_id%%molecule_atoms )%in% edge_atoms
                        ,0
                        ,floor(atom_id/molecule_atoms)*molecule_monomers
                        +ceiling((atom_id%%molecule_atoms+ifelse(atom_id%%molecule_atoms>16,-1,0))/monomer_atoms)
      )
      return(monomer.id)
    }
    
  }
  ##########
  tot_atom_num=atom.300.1_fread[,.(atom.id)]%>%max()
  
  timestep=dim(atom.300.1_fread)[1]/tot_atom_num
  
  paste("The number of timestep:",timestep,sep="")%>%message()
  
  ################################################################################################
  
  # calculate the xu,yu,zu unscaled position 
  
  ## first give the bounds of box on the x,y,z direction. 
  
  
  
  xlength=xhigh-xlow
  ylength=yhigh-ylow
  zlength=zhigh-zlow
  
  
  atom.300.1_fread[,`:=`(xu=xs*xlength+xlow+ix*xlength,yu=ys*ylength+ylow+iy*ylength,zu=zs*zlength+zlow+iz*zlength)]
  
  #####################################################################
  
  # add mol variable to original data frame
  
  timeRecordB()
  atom.300.1_fread[,mol:=mol_id_gen(atom_id = atom.id,molecule_atoms=molecule_atoms)]
  timeRecordB(output_message = "add mol variable data.table")
  ################################################################################################
  # add time variable for easier data handling and safe guarding any unwanted ordering to change the ordering of time in the data. 
  timeRecordB()
  atom.300.1_fread[,time_step:=seq(1,timestep,by=1)%>%rep(times=tot_atom_num)%>%sort()]
  timeRecordB(output_message = "add time variable data.table")
  
  # add mass variable to original data frame  
  timeRecordB()
  atom.300.1_fread[,mass:=decode(x = type, search = atom_type
                                 , replace = atom_type_mass)]
  timeRecordB(output_message = "add mass variable data.table")
  gc()
  
  
  # add monomer id 
  timeRecordB()
  atom.300.1_fread[,monomer.id:=monomer_gen(atom_id=atom.id,molecule_atoms=molecule_atoms,molecule_monomers=molecule_monomers,monomer_atoms=monomer_atoms)]
  timeRecordB(output_message = "add monomer id variable data.table")
  gc()
  
  # mass of monomer 
  monomer_mass=atom.300.1_fread[time_step==1,]%>%.[monomer.id==1,.(mass)]%>%sum
  
  
  MSD_g1_matrix=function(data,timestep,num_monomer,monomer_mass){
    MSD_g1_empty_matrix=matrix(NA,ncol=timestep,nrow=num_monomer)
    ###############################################################
    # calculate the center of mass for every timestep and every monomer 
    center_of_mass=data[,.(xu=sum(xu*mass)/monomer_mass,yu=sum(yu*mass)/monomer_mass,zu=sum(zu*mass)/monomer_mass),by=.(monomer.id,time_step)]
    
    ###############################################################
    for(j in 1:num_monomer){
      MSD_g1_empty_matrix[j,]=center_of_mass[monomer.id==j,]%>%MSD
    }
    return(MSD_g1_empty_matrix)
  } 
  
  timeRecordB()
  
  MSD.matrix=MSD_g1_matrix(atom.300.1_fread,timestep=timestep,num_monomer= num_mol*molecule_monomers,monomer_mass=monomer_mass)
  
  timeRecordB(output_message = "MSD for center of mass")
  
  gc()
  
  # Calculate the averaged MSD over all molecules. 
  
  MSD.matrix%>%colMeans()%>%write.table(file="MSD.g1.colmean.1.txt", sep=",")
  
  
  
  ## Add NGP calculation here. 
  timeRecordB()
  
  NGP.COM=(0.6)*colMeans(MSD.matrix^2)/(colMeans(MSD.matrix))^2-1
  
  NGP.COM%>%write.table(file="NGP.g1.1.txt", sep=",")
  
  timeRecordB(output_message = "NGP for center of mass of monomer")
  
  
  return( timeRecordR(ignore=0.1)%>%filter(output_message!="None")%>%select(output_message,run_time) )
}


# example use # put the timeRecord function loading outside the function in case internet conenction is down during the calculation. 
# source("https://raw.githubusercontent.com/edwardcooper/mlmodel_select/master/timeRecord_functions.R")

# MSD_scaled_g1_one_temp(path="~/Dropbox/lammps/PMMA_long/atom300",filename="atom.300_long2",polymer="PMMA_long"
#                        ,num_mol=64,molecule_atoms=602,molecule_monomers=40 ,monomer_atoms=15
#                        ,atom_type=1:10,atom_type_mass=c(1.0079,12.011,12.011,12.011,15.9999,15.9999,12.011,12.011,1.0079,12.011)
#                        ,xhigh=72.3695
#                        ,xlow=2.86271
#                        ,yhigh=75.7675
#                        ,ylow=-0.535315
#                        ,zhigh=73.6408
#                        ,zlow=1.59139)

## echo the current calculation and percentage of entire calculation.  
MSD_scaled_g1=function(Path="~/Dropbox/lammps/",polymer="PMMA_long",temperatures=seq(300,620,by=20)
                       ,num_mol=64,molecule_atoms=602,molecule_monomers=40 ,monomer_atoms=15
                       ,atom_type=1:10,atom_type_mass=c(1.0079,12.011,12.011,12.011,15.9999,15.9999,12.011,12.011,1.0079,12.011)
                       ,xhigh=72.3695
                       ,xlow=2.86271
                       ,yhigh=75.7675
                       ,ylow=-0.535315
                       ,zhigh=73.6408
                       ,zlow=1.59139
){
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
    filename=paste("atom.",temperatures[i],"_long2",sep="")
    
    # calculation for MSD
    MSD_scaled_g1_one_temp(path=path,filename =filename,num_mol=num_mol,molecule_atoms=molecule_atoms,molecule_monomers=molecule_monomers,monomer_atoms=monomer_atoms
                           ,atom_type=atom_type,atom_type_mass=atom_type_mass
                           ,xhigh=xhigh,xlow=xlow,yhigh=yhigh,ylow=ylow,zhigh=zhigh,zlow=zlow )
    
    # echo end of calculation
    paste("End calculation of temperature:",temperatures[i],sep="")%>%message
    paste(i,"/",length(temperatures))%>%message
    gc()
  }
  
  return( timeRecordR(ignore=0.1)%>%filter(output_message!="None")%>%select(output_message,run_time) )
}

MSD_scaled_g1(Path="~/Dropbox/lammps/",polymer="PMMA_long",temperatures=seq(300,620,by=20)
              ,num_mol=64,molecule_atoms=602,molecule_monomers=40 ,monomer_atoms=15
              ,atom_type=1:10,atom_type_mass=c(1.0079,12.011,12.011,12.011,15.9999,15.9999,12.011,12.011,1.0079,12.011)
              ,xhigh=72.3695
              ,xlow=2.86271
              ,yhigh=75.7675
              ,ylow=-0.535315
              ,zhigh=73.6408
              ,zlow=1.59139)
