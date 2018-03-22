######################################################

# short time results 

# g2 Tmean for scaled PS_20 
source('~/Dropbox/lammps/MSD_scaled_g2_Tmean.R', echo=TRUE)


MSD_scaled_g2_Tmean(Path="~/Dropbox/lammps/",polymer="PS_20",temperatures=seq(200,520,by=20)
                   ,molecule_atoms=645,num_mol=40,atom_type=1:6,atom_type_mass=c(12.011,1.0079,12.011,12.011,12.011,1.0079)
                   ,core_num=4)
# Done 


#################### 


# g1 T mean for unscaled PMMA_big 
source('~/Dropbox/lammps/MSD_unscaled_g1_Tmean.R', echo=TRUE)


MSD_unscaled_g1_Tmean(Path="~/Dropbox/lammps/",polymer="PMMA_big",temperatures=seq(300,600,by=20)
                      ,num_mol=64
                      ,molecule_atoms=602
                      ,molecule_monomers=40
                      ,monomer_atoms=15
                      ,atom_type=1:10
                      ,atom_type_mass=c(1.0079,12.011,12.011,12.011,15.9999,15.9999,12.011,12.011,1.0079,12.011)
                      ,core_num=4)
 

# g2 Tmean for unscaled PMMA_big
source('~/Dropbox/lammps/MSD_unscaled_g2_Tmean.R', echo=TRUE)

MSD_unscaled_g2_Tmean(Path="~/Dropbox/lammps/",polymer="PMMA_big",temperatures=seq(300,600,by=20)
                               ,molecule_atoms=602,num_mol=64,atom_type=1:10
                               ,atom_type_mass=c(1.0079,12.011,12.011,12.011,15.9999,15.9999,12.011,12.011,1.0079,12.011)
                               ,core_num=4)



######################################################

# long time results 


## PS_20_long
# g0 T mean for scaled PS_20_long
source('~/Dropbox/lammps/MSD_scaled_g0_Tmean.R',echo=TRUE)
MSD_scaled_g0_Tmean(Path="~/Dropbox/lammps/",polymer="PS_20_long",temperatures=seq(200,520,by=20), core_num = 10)

# g1 T mean for scaled PS_20_long
source('~/Dropbox/lammps/MSD_scaled_g1_Tmean.R', echo=TRUE)

MSD_scaled_g1_Tmean(Path="~/Dropbox/lammps/",polymer="PS_20_long",temperatures=seq(200,520,by=20)
                    ,num_mol=40
                    ,molecule_atoms=645
                    ,molecule_monomers=40
                    ,monomer_atoms=16
                    ,atom_type=1:6
                    ,atom_type_mass=c(12.011,1.0079,12.011,12.011,12.011,1.0079)
                    ,core_num=10)

# g2 T mean for scaled PS_20_long

source('~/Dropbox/lammps/MSD_scaled_g2_Tmean.R', echo=TRUE)

MSD_scaled_g2_Tmean(Path="~/Dropbox/lammps/",polymer="PS_20_long",temperatures=seq(200,520,by=20)
                    ,molecule_atoms=645,num_mol=40,atom_type=1:6,atom_type_mass=c(12.011,1.0079,12.011,12.011,12.011,1.0079)
                    ,core_num=10)

#################### 

## PMMA_long

# g0 T mean for scaled PMMA_long
source('~/Dropbox/lammps/MSD_scaled_g0_Tmean.R',echo=TRUE)
MSD_scaled_g0_Tmean(Path="~/Dropbox/lammps/",polymer="PMMA_long",temperatures=seq(300,640,by=20), core_num = 10)

# g1 T mean for scaled PMMA_long 

source('~/Dropbox/lammps/MSD_scaled_g1_Tmean.R', echo=TRUE)

MSD_scaled_g1_Tmean(Path="~/Dropbox/lammps/",polymer="PMMA_long",temperatures=seq(300,640,by=20)
                      ,num_mol=64
                      ,molecule_atoms=602
                      ,molecule_monomers=40
                      ,monomer_atoms=15
                      ,atom_type=1:10
                      ,atom_type_mass=c(1.0079,12.011,12.011,12.011,15.9999,15.9999,12.011,12.011,1.0079,12.011)
                      ,core_num=10)


# g2 Tmean for unscaled PMMA_long
source('~/Dropbox/lammps/MSD_scaled_g2_Tmean.R', echo=TRUE)

MSD_scaled_g2_Tmean(Path="~/Dropbox/lammps/",polymer="PMMA_long",temperatures=seq(300,640,by=20)
                      ,molecule_atoms=602,num_mol=64,atom_type=1:10
                      ,atom_type_mass=c(1.0079,12.011,12.011,12.011,15.9999,15.9999,12.011,12.011,1.0079,12.011)
                      ,core_num=10)


