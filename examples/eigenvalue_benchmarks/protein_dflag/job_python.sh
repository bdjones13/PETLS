#!/bin/bash --login

#SBATCH --time=1:00:00
#SBATCH --nodes=1
## SBATCH --ntasks=1
#SBATCH -c 1
#SBATCH --mem=32g
#SBATCH --job-name=scipy
#SBATCH --output=/mnt/home/jones657/Documents/PersistentLaplacians/examples/eigenvalue_benchmarks/protein_dflag/outfiles/%x-%A-%a.out
                 
#SBATCH --array=0-99

source /mnt/home/jones657/Documents/venv_petls/bin/activate
cd /mnt/home/jones657/Documents/PersistentLaplacians/examples/eigenvalue_benchmarks/protein_dflag/


ALGS=(
    "scipy.sparse.linalg.eigsh.largest"
    "scipy.sparse.linalg.eigsh.smallest.shifted"
    "scipy.sparse.linalg.eigsh.smallest.shifted.tol0"
    "scipy.linalg.eigvalsh"
)

for alg in "${ALGS[@]}"; do
    echo "Running $alg"
    python3 compute_eigenvalues.py $alg $SLURM_ARRAY_TASK_ID
done
