#!/bin/bash  
#SBATCH --partition=SP2
#SBATCH --ntasks=10 		# number of tasks / mpi processes
#SBATCH --cpus-per-task=1 	# Number OpenMP Threads per process
#SBATCH -J cobaya_test1
#SBATCH --time=8-00:00:00

#E-mail notifications:
#SBATCH --mail-type=ALL
#SBATCH --mail-user=meirellesarthur@usp.br

#OpenMP settings:
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

echo $SLURM_JOB_ID		#ID of job allocation
echo $SLURM_SUBMIT_DIR		#Directory job where was submitted
echo $SLURM_JOB_NODELIST	#File containing allocated hostnames
echo $SLURM_NTASKS		#Total number of cores for job

module swap openmpi mpich
module swap gnu gnu7

cd /temporario2/10300487/MscProject/cross-correlation/theoretical/cobaya/

mpirun -n 2 python pyctg_likelihood.py > scriptrun_test/out.log 2>&1

