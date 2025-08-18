#!/bin/bash --login

#SBATCH --time=0:05:00
#SBATCH --nodes=1
## SBATCH --ntasks=1
#SBATCH -c 1
#SBATCH --mem=16g
#SBATCH --job-name=create_toruses
#SBATCH --output=/mnt/home/jones657/Documents/PersistentLaplacians/examples/eigenvalue_benchmarks/outfiles/%x-%A-%a.out
                 
#SBATCH --array=0-99

source /mnt/home/jones657/Documents/venv_petls/bin/activate
cd /mnt/home/jones657/Documents/PersistentLaplacians/examples/eigenvalue_benchmarks/
python3 create_toruses.py $SLURM_ARRAY_TASK_ID
