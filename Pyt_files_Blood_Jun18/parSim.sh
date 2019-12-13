#!/bin/bash -l

#SBATCH
#SBATCH --job-name=MyJob
#SBATCH --time=3:0:0
#SBATCH --partition=parallel
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=24
# number of cpus (threads) per task (process)
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=end
#SBATCH --mail-user=klahoue1@jhu.edu

module load python/3.6-anaconda

python Corse_tau_leap_less_param_Blood > $SLURM_JOBID.txt

