#!/bin/bash -l

#SBATCH
#SBATCH --job-name=MyJob
#SBATCH --time=2:0:0
#SBATCH --partition=parallel
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=24
# number of cpus (threads) per task (process)
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=end
#SBATCH --mail-user=klahoue1@jhu.edu

module load python/3.6-anaconda

python Coarse_tau_leap_less_paramNone > $SLURM_JOBID.txt

