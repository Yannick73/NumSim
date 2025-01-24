#!/bin/bash

#

#SBATCH --job-name=NumSim1N

#SBATCH --output=profiling1N.txt

#

#SBATCH --ntasks=48

#SBATCH --ntasks-per-node=48

#SBATCH --time=12:00:00

#SBATCH --exclusive


module use /usr/local.nfs/sgs/modulefiles

module load gcc/10.2

module load openmpi/3.1.6-gcc-10.2

module load vtk/9.0.1

module load cmake/3.18.2

echo "NumSim Ex2 profiling with Scenario3-6"

for i in 48 32 20 16 12 8 6 4 2 1
do
	echo "\nSIMULATE WITH $i RANK\n"
	srun -n $i ./build/numsim_parallel input/scenario3_6.txt
done
