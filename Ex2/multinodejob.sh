#!/bin/bash

#

#SBATCH --job-name=NumSim2N

#SBATCH --output=profiling2N.txt

#

#SBATCH --ntasks=96

#SBATCH --ntasks-per-node=48

#SBATCH --time=12:00:00

#SBATCH --exclusive


module use /usr/local.nfs/sgs/modulefiles

module load gcc/10.2

module load openmpi/3.1.6-gcc-10.2

module load vtk/9.0.1

module load cmake/3.18.2

echo "NumSim Ex2 profiling with Scenario3-6"

for i in 48 64 80 96
do
	echo "\nSIMULATE WITH $i RANK\n"
	srun -n $i ./build/numsim_parallel input/scenario3_6.txt
done
