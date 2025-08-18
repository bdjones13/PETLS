#!/bin/bash --login

#SBATCH --time=1:00:00
#SBATCH --nodes=1
## SBATCH --ntasks=1
#SBATCH -c 1
#SBATCH --mem=32g
#SBATCH --job-name=spectra_smallest
#SBATCH --output=/mnt/home/jones657/Documents/PersistentLaplacians/examples/eigenvalue_benchmarks/protein_dflag/outfiles/%x-%A-%a.out
                 
#SBATCH --array=0-99


cd /mnt/home/jones657/Documents/PersistentLaplacians/examples/eigenvalue_benchmarks/protein_dflag/

# ./build/cpp_eigs eigen.selfadjoint  $SLURM_ARRAY_TASK_ID
# echo "done with eigen.selfadjoint"

# ./build/cpp_eigs eigen.eigensolver  $SLURM_ARRAY_TASK_ID
# echo "done with eigen.eigensolver"

# ./build/cpp_eigs spectra.dense.largest  $SLURM_ARRAY_TASK_ID
# echo "done with spectra.dense.largest"

./build/cpp_eigs spectra.inverse.smallest  $SLURM_ARRAY_TASK_ID
echo "done with spectra.inverse.smallest"