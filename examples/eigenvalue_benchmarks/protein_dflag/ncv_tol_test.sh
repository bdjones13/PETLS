#!/bin/bash --login

#SBATCH --time=12:00:00
#SBATCH --nodes=1
## SBATCH --ntasks=1
#SBATCH -c 1
#SBATCH --mem=32g
#SBATCH --job-name=tol13
#SBATCH --output=/mnt/home/jones657/Documents/PersistentLaplacians/examples/eigenvalue_benchmarks/protein_dflag/outfiles/%x-%A-%a.out
                 
## SBATCH --array=0-99

source /mnt/home/jones657/Documents/venv_petls/bin/activate
cd /mnt/home/jones657/Documents/PersistentLaplacians/examples/eigenvalue_benchmarks/protein_dflag


# python3 ncv_tol_test.py 0
# python3 ncv_tol_test.py 1e-4
# python3 ncv_tol_test.py 1e-6
# python3 ncv_tol_test.py 1e-11
# python3 ncv_tol_test.py 1e-12
python3 ncv_tol_test.py 1e-13

# python3 ncv_tol_test.py 1e-14
# python3 ncv_tol_test.py 1e-10