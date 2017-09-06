#!/bin/bash
#SBATCH -N 2 -n 12
#SBATCH --time=100:00:00

module load lammps

cp input_script data.in.txt  $PFSDIR
cd $PFSDIR
mpiexec lmp2015 < input_script
cp log.lammps $SLURM_SUBMIT_DIR
