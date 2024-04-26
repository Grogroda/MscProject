#!/bin/bash

#SBATCH --partition=SP2
#SBATCH --ntasks=1 		# number of tasks / mpi processes
#SBATCH --cpus-per-task=9 	# Number OpenMP Threads per process
#SBATCH -J band1_profile 
#SBATCH --time=8-00:00:00

#E-mail notifications:
#SBATCH --mail-type=ALL
#SBATCH --mail-user=meirellesarthur@usp.br

#OpenMP settings:
#export OMP_NUM_THREADS=1
#export MKL_NUM_THREADS=1
#export OMP_PLACES=threads
#export OMP_PROC_BIND=spread

echo $SLURM_JOB_ID		#ID of job allocation
echo $SLURM_SUBMIT_DIR		#Directory job where was submitted
echo $SLURM_JOB_NODELIST	#File containing allocated hostnames
echo $SLURM_NTASKS		#Total number of cores for job

module swap openmpi mpich
module swap gnu gnu7

source /temporario2/10300487/miniconda3/bin/activate

#cd /temporario2/10300487/MscProject/cross-correlation/theoretical/src/

#./compile_lib.sh

#echo "Lib compiled!"

cd /temporario2/10300487/MscProject/cross-correlation/theoretical/cobaya/

python3 profile_OmegaM.py -n 9 -b 1 -N 20000000 -c 'ctg' --npoints 50 --Omin 0.1 --Omax 0.5 -i 1 > like_profiles/band1_profile_parallel1.log 2>&1 

#python3 profile_OmegaM.py -h -b -n -N -c --npoints --Omin --Omax
#h=help, b=band (1,2,3,4 or 'min'), n=number of processes
#N=number of MC calls, -c='ctg' or 'cgg', npoints=number of OmegaM's
#Omin=Omega_min and Omax=Omega_max
