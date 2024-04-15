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

echo "Lib compiled!"

cd /temporario2/10300487/MscProject/cross-correlation/theoretical/cobaya/

python3 profile_OmegaM.py -n 9 -b 1 > like_profiles/band1_profile_parallel.log 2>&1 #number after script name=2MASS band
