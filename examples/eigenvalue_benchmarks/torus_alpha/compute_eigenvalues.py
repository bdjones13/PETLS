# Note: run create_toruses.ipynb first to generate the required L matrices in the torus_L directory

from io import StringIO
from scipy.io import mmread
import sys
import os
import pandas as pd

import timeit


from scipy.linalg import eig, eigh
import numpy as np
import scipy.sparse.linalg

def get_filtrations(delta):
    start = 0.0
    end = 5.0
    filtrations = [start+i*delta for i in range(int((end-start)/delta) + 1)]
    return filtrations


def sparse_wrapper(L,alg="eigs",which="LM",sigma=None):

    k=10
    simplices = L.shape[0]
    k = min(k, simplices - 1) if simplices > 1 else 1
    
    
    
    ncv = min(simplices, max(30, 3*k))

    if alg == "eigs":
        try:
            return np.array(sorted(scipy.sparse.linalg.eigs(L, k=k, which=which, ncv=ncv, sigma=sigma)[0].real))
        except scipy.sparse.linalg.ArpackNoConvergence as e:
            print("arpack non-convergence error")
            return np.array(sorted(scipy.linalg.eigvals(L.todense()).real))
    elif alg == "eigsh":
        try:
            return scipy.sparse.linalg.eigsh(L, k=k, which=which, ncv=ncv, sigma=sigma)
        except scipy.sparse.linalg.ArpackNoConvergence as e:
            print("arpack non-convergence error")
            return np.array(sorted(scipy.linalg.eigvalsh(L.todense())))




def do_test(algorithm, replicate):
    filtrations = get_filtrations(delta=0.1)
    dims = [0, 1, 2]
    
    
    times = []
    # Loop through each filtration and dimension
    for filtration in filtrations:
            for dim in dims:
                    # read file
                    filename = f"./torus_L/dim{dim}_a{filtration:.2f}_b{filtration+0.1:.2f}_r{replicate}.mkt"
                    L = mmread(filename)
                    # print(f"File {filename} has shape {L.shape}")
                    simplices = L.shape[0]
                    if simplices == 0:
                         times.append([dim, filtration, filtration+0.1, simplices, 0.0])
                         continue
                    
                    if "sparse" in algorithm:

                        if simplices == 1:
                            times.append([dim, filtration, filtration+0.1, simplices, -1])
                            continue
                        
                        if algorithm == "scipy.sparse.linalg.eigs":
                            t = timeit.timeit(lambda: sparse_wrapper(L), number=1)

                        elif algorithm == "scipy.sparse.linalg.eigsh.largest":
                            t = timeit.timeit(lambda: sparse_wrapper(L, alg="eigsh"), number=1)

                        elif algorithm == "scipy.sparse.linalg.eigs.smallest":
                            t = timeit.timeit(lambda: sparse_wrapper(L, which="SM"), number=1)

                        elif algorithm == "scipy.sparse.linalg.eigsh.smallest":
                            t = timeit.timeit(lambda: sparse_wrapper(L, alg="eigsh",which="SM"), number=1)

                        elif algorithm == "scipy.sparse.linalg.eigs.smallest.shifted":
                            t = timeit.timeit(lambda: sparse_wrapper(L, sigma=1e-3), number=1)

                        elif algorithm == "scipy.sparse.linalg.eigsh.smallest.shifted":
                            t = timeit.timeit(lambda: sparse_wrapper(L,alg="eigsh", sigma=1e-3), number=1)

                    else:    
                        L = L.todense()
                        #compute eigenvalues
                        if algorithm == "scipy.linalg.eigvals":
                            t = timeit.timeit(lambda: eig(L), number=1)
                        elif algorithm == "scipy.linalg.eigvalsh":
                            t = timeit.timeit(lambda: eigh(L, eigvals_only=True), number=1)
                        # elif algorithm == "numpy.linalg.eigvals":
                        #     t = timeit.timeit(lambda: np.linalg.eigvals(L), number=1)
                        # elif algorithm == "numpy.linalg.eigvalsh":
                        #     t = timeit.timeit(lambda: np.linalg.eigvalsh(L), number=1)
                        
                    print(f"Time taken for {algorithm}: {t:.4f} seconds",flush=True)
        
                        # add to list
                    times.append([dim,filtration, filtration+0.1, simplices, t])
    
    # output results to csv
    df = pd.DataFrame(times, columns=["dim", "a", "b", "simplices",algorithm])
    foldername = f"profiles/{algorithm}"
    if not os.path.exists(foldername):
        os.mkdir(foldername)
    fname = f"{foldername}/r{replicate}.csv"
    df.to_csv(fname, index=False, header=True, columns=["dim","a","b","simplices",f"{algorithm}"])
    print(f"Results saved to {fname}", flush=True)

if __name__ == "__main__":
    algorithm=sys.argv[1]
    replicate=sys.argv[2]
    # replicate = 0
    # algorithm = "scipy.linalg.eigvalsh"
    do_test(algorithm, replicate)    