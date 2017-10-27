library(magrittr)
# the last atom is defined as 0. 
monomer_gen_PS=function(atom_id,molecule_atoms=645,molecule_monomers=40,monomer_atoms=16,edge_atoms=c(17,642,643,644,0)){

  monomer.id=ifelse((atom_id%%molecule_atoms )%in% edge_atoms
                    ,0
                    ,floor(atom_id/molecule_atoms)*molecule_monomers
                    +ceiling((atom_id%%molecule_atoms+ifelse(atom_id%%molecule_atoms>16,-1,0))/monomer_atoms)
                    )
  return(monomer.id)
}

data.frame(id=monomer_gen_PS(atom_id =1:34 ), atom.id=1:34  )

data.frame(id=monomer_gen_PS(atom_id =610:645), atom.id=610:645 )

monomer_gen_PS(atom_id =1:25800)%>%table()%>%table()




## Monomer id gen function test for PMMA_big

# the last atom is defined as 0. 
monomer_gen=function(atom_id,molecule_atoms=602,molecule_monomers=40,monomer_atoms=15){
  # if it is on the either end of polymer, the monomer is defined as monomer 0. 
  # if the atom is not on the end, then first calculate the molecule number and multiply it by 40 since there is 40 monomers in each molecule.
  # then add the the monomer number it has in this molecule. Need to deduct the first atom off the molecule, divided by the number of atoms in a monomer and get a integer
  # then you have the monomer id. 
  monomer.id=ifelse(atom_id%%molecule_atoms %in% c(0,1)
                    ,0
                    ,floor(atom_id/molecule_atoms)*molecule_monomers+ceiling((atom_id%%molecule_atoms-1)/monomer_atoms) 
                    )
  return(monomer.id)
}

data.frame(id=monomer_gen(atom_id =1:34 ), atom.id=1:34  )

data.frame(id=monomer_gen(atom_id =573:602 ), atom.id=573:602  )

monomer_gen(atom_id =1:(602*64))%>%table%>%table()
           