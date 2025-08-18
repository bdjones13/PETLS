# Note: run create_toruses.ipynb first to generate the required L matrices in the torus_L directory

from io import StringIO
from scipy.io import mmread
import sys
import os
import pandas as pd

import timeit
import time

from scipy.linalg import eigvalsh
import numpy as np
import scipy.sparse.linalg


def get_filtrations(delta):
    start = 0.0
    end = 10
    filtrations = [start+i*delta for i in range(int((end-start)/delta) + 1)]
    return filtrations


def sparse_wrapper(L,alg="eigsh",which="LM",sigma=None,tol=0):

    k=10
    simplices = L.shape[0]
    
    ncv = min(simplices, 60)

    if alg == "eigsh":
        try:
            k = min(k, simplices - 1) if simplices > 1 else 1
            return scipy.sparse.linalg.eigsh(L, k=k, which=which, ncv=ncv, sigma=sigma, tol=tol, maxiter=2000)
        except scipy.sparse.linalg.ArpackNoConvergence as e:
            # print("arpack non-convergence error")
            eigs = np.array(sorted(scipy.linalg.eigvalsh(L.todense())))
            raise e




def do_test(algorithm, replicate):
    filtrations = get_filtrations(delta=0.1)
    # dims = [0, 1, 2]
    dims = [1, 2]
    
    times = []
    # Loop through each filtration and dimension
    for filtration in filtrations:
            for dim in dims:
                    failed = False
                    # read file
                    filename = f"./protein_L/dim{dim}_a{filtration:.2f}_b{filtration+0.1:.2f}.mkt"
                    # print(f"Try to read file: {filename}")
                    L = mmread(filename)
                    # print(f"File {filename} has shape {L.shape}")
                    simplices = L.shape[0]
                    if simplices == 0 or L.count_nonzero() == 0:
                         times.append([dim, filtration, filtration+0.1, simplices, 0.0])
                         continue
                    
                    if "sparse" in algorithm:
                        backup_timer = time.time()
                        if simplices == 1:
                            times.append([dim, filtration, filtration+0.1, simplices, -1])
                            continue
                        try:
                            if algorithm == "scipy.sparse.linalg.eigsh.largest":
                                t = timeit.timeit(lambda: sparse_wrapper(L, alg="eigsh",tol=0), number=1)
                            elif algorithm == "scipy.sparse.linalg.eigsh.smallest.shifted":
                                t = timeit.timeit(lambda: sparse_wrapper(L,alg="eigsh", sigma=1e-2,tol=1e-6), number=1)
                            elif algorithm == "scipy.sparse.linalg.eigsh.smallest.shifted.tol0":
                                t = timeit.timeit(lambda: sparse_wrapper(L,alg="eigsh", sigma=1e-2,tol=0.0), number=1)
                        except Exception as e:
                            print(f"Exception: {e}")
                            failed = True
                            t = time.time()-backup_timer
                            # raise e
                    else:    
                        L = L.todense()

                        if algorithm == "scipy.linalg.eigvalsh":
                            t = timeit.timeit(lambda: scipy.linalg.eigvalsh(L), number=1)

                    print(f"Time taken for {algorithm} at d={dim}, filtration={filtration}: {t:.4f} seconds",flush=True)
        
                        # add to list
                    times.append([dim,filtration, filtration+0.1, simplices, t, failed])
                    failed = False
    
    # output results to csv
    df = pd.DataFrame(times, columns=["dim", "a", "b", "simplices",algorithm,f"{algorithm}_failed"])
    foldername = f"scratch_profiles/{algorithm}"
    if not os.path.exists(foldername):
        os.mkdir(foldername)
    fname = f"{foldername}/r{replicate}.csv"
    df.to_csv(fname, index=False, header=True, columns=["dim","a","b","simplices",f"{algorithm}",f"{algorithm}_failed"])
    print(f"Results saved to {fname}", flush=True)

if __name__ == "__main__":
    algorithm=sys.argv[1]
    replicate=sys.argv[2]
    # if len(sys.argv) > 3:
    #     tol = float(sys.argv[3])
    # else:
    #     tol = 0
    # replicate = 0
    # algorithm = "scipy.linalg.eigvalsh"
    do_test(algorithm, replicate)#, tol=tol)    