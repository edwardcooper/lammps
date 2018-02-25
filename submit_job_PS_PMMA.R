
source("https://raw.githubusercontent.com/edwardcooper/lammps/master/MSD_g0.R")
source("https://raw.githubusercontent.com/edwardcooper/lammps/master/MSD_g1.R")
source("https://raw.githubusercontent.com/edwardcooper/lammps/master/MSD_g2.R")





MSD_g0(Path="~/Dropbox/lammps/",polymer="PMMA_big",temperatures=seq(300,600,by=20),timestep=5001)

MSD_g1(Path="~/Dropbox/lammps/",polymer="PMMA_big",temperatures=seq(300,600,by=20),timestep=5001
          ,num_mol=64
          ,molecule_atoms=602,molecule_monomers=40,monomer_atoms=15
          ,atom_type=1:10,atom_type_mass=c(1.0079,12.011,12.011,12.011,15.9999,15.9999,12.011,12.011,1.0079,12.011))

MSD_g2(Path="~/Dropbox/lammps/",polymer="PMMA_big",temperatures=seq(300,600,by=20),timestep=5001
       ,num_mol=64,atom_type=1:10,atom_type_mass=c(1.0079,12.011,12.011,12.011,15.9999,15.9999,12.011,12.011,1.0079,12.011))




MSD_g0(Path="~/Dropbox/lammps/",polymer="PS",temperatures=seq(200,600,by=50),timestep=5001)

MSD_g1(Path="~/Dropbox/lammps/",polymer="PS",temperatures=seq(200,600,by=50),timestep=5001
      ,num_mol=40
      ,molecule_atoms=645,molecule_monomers=40,monomer_atoms=16
      ,atom_type=1:6,atom_type_mass=c(12.011,1.0079,12.011,12.011,12.011,1.0079))

MSD_g2(Path="~/Dropbox/lammps/",polymer="PS",temperatures=seq(200,600,by=50),timestep=5001
         ,num_mol=40,atom_type=1:6,atom_type_mass=c(12.011,1.0079,12.011,12.011,12.011,1.0079))



## g0 MSD: you could use the same function for g0 calculation since you just calculate the MSD of each atom and average them all. 

## g1 MSD: you need different functions to add monomer id (edge atoms are added as monomer.id =0), but center of mass calculation and MSD calculation is the same once you have the monomer id.

## g2 MSD: you need the same function but different atom_type and masses to do it. (plus, the problem with mol starts at 0 or 1)/Solved




# test each function run before you launch to do it for all temperatures. 

source("https://raw.githubusercontent.com/edwardcooper/lammps/master/MSD_scaled_g0_revised.R")
source("https://raw.githubusercontent.com/edwardcooper/lammps/master/MSD_scaled_g1_revised.R")
source("https://raw.githubusercontent.com/edwardcooper/lammps/master/MSD_scaled_g2_revised.R")

MSD_scaled_g0(Path="~/Dropbox/lammps/",polymer="PS_20",temperatures=seq(200,520,by=20)
              )
MSD_scaled_g1(Path="~/Dropbox/lammps/",polymer="PS_20",temperatures=seq(200,520,by=20)
                     ,molecule_atoms=645,num_mol=40,molecule_monomers=40
                     ,monomer_atoms=16,atom_type=1:6,atom_type_mass=c(12.011,1.0079,12.011,12.011,12.011,1.0079) )
MSD_scaled_g2(Path="~/Dropbox/lammps/",polymer="PS_20",temperatures=seq(200,520,by=20)
             ,molecule_atoms=645,num_mol=40,atom_type=1:6,atom_type_mass=c(12.011,1.0079,12.011,12.011,12.011,1.0079))
                   