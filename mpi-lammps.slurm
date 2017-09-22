#!/bin/bash
#SBATCH -N 10 -n 60
#SBATCH --time=100:00:00
#SBATCH --mail-user=zxl220@case.edu
module load lammps

cp input_script data.in.txt  $PFSDIR
cd $PFSDIR
mpiexec lmp2015 < input_script
cp log.lammps $SLURM_SUBMIT_DIR
