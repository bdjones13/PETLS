#!/bin/bash --login

#SBATCH --time=1:00:00
#SBATCH --nodes=1
## SBATCH --ntasks=1
#SBATCH -c 1
#SBATCH --mem=32g
#SBATCH --job-name=large_alpha
#SBATCH --output=/mnt/home/jones657/Documents/PersistentLaplacians/examples/large_complex/outfiles/%x-%A-%a.out
                 
##SBATCH --array=

source /mnt/home/jones657/Documents/venv_petls/bin/activate
cd /mnt/home/jones657/Documents/PersistentLaplacians/examples/large_complex

python3 large_alpha.py
