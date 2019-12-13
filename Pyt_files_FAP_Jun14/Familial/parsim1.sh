#!/bin/bash -l

#SBATCH
#SBATCH --job-name=MyJob
#SBATCH --time=24:0:0
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=24
# number of cpus (threads) per task (process)
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=end
#SBATCH --mail-user=klahoue1@jhu.edu

module load python/2.7.10

python Tau_leap_Colon_Familial.py > $SLURM_JOBID.txt
