#!/bin/bash --login

#SBATCH --time=12:00:00
#SBATCH --nodes=1
## SBATCH --ntasks=1
#SBATCH -c 1
#SBATCH --mem=32g
#SBATCH --job-name=python_eigs
#SBATCH --output=/mnt/home/jones657/Documents/PersistentLaplacians/examples/eigenvalue_benchmarks/outfiles/%x-%A-%a.out
                 
#SBATCH --array=63 #0-99

source /mnt/home/jones657/Documents/venv_petls/bin/activate
cd /mnt/home/jones657/Documents/PersistentLaplacians/examples/eigenvalue_benchmarks/


ALGS=(
    "scipy.linalg.eigvalsh"
    "scipy.linalg.eigvals"
    "scipy.sparse.linalg.eigs"
    "scipy.sparse.linalg.eigsh"    
    "scipy.sparse.linalg.eigs.smallest"
    "scipy.sparse.linalg.eigsh.smallest"
    "scipy.sparse.linalg.eigs.smallest.shifted"
    "scipy.sparse.linalg.eigsh.smallest.shifted"
)

for alg in "${ALGS[@]}"; do
    echo "Running $alg"
    python3 compute_eigenvalues.py $alg $SLURM_ARRAY_TASK_ID
done
