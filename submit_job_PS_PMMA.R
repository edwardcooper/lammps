
source("https://raw.githubusercontent.com/edwardcooper/lammps/master/MSD_g0.R")
source("https://raw.githubusercontent.com/edwardcooper/lammps/master/MSD_g1.R")
source("https://raw.githubusercontent.com/edwardcooper/lammps/master/MSD_g2.R")





MSD_g0(Path="~/Dropbox/lammps/",polymer="PMMA_big",temperatures=seq(300,600,by=20),timestep=5001)
MSD_g1(Path="~/Dropbox/lammps/",polymer="PMMA_big",temperatures=seq(300,600,by=20),timestep=5001)
MSD_g2(Path="~/Dropbox/lammps/",polymer="PMMA_big",temperatures=seq(300,600,by=20),timestep=5001)


MSD_g0(Path="~/Dropbox/lammps/",polymer="PS",temperatures=seq(200,600,by=50),timestep=5001)





# Need to modify the MSD_g1() functions to deal with PS monomer MSD.


# Need some more data on atom_type and mass + number of molecules in total+ adjust the starting mol number. 
MSD_g2(Path="~/Dropbox/lammps/",polymer="PS",temperatures=seq(200,600,by=50),timestep=5001
         ,num_mol=40,atom_type=1:6,atom_type_mass=c(12.011,1.0079,12.011,12.011,12.011,1.0079))



## g0 MSD: you could use the same function for g0 calculation since you just calculate the MSD of each atom and average them all. 

## g1 MSD: you need different functions to add monomer id (edge atoms are added as monomer.id =0), but center of mass calculation and MSD calculation is the same once you have the monomer id.

## g2 MSD: you need the same function but different atom_type and masses to do it. (plus, the problem with mol starts at 0 or 1)/Solved