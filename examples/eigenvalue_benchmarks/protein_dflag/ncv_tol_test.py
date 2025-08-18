import scipy
import numpy as np
from scipy.io import mmread
import time
import timeit
import pandas as pd
import sys

def get_filtrations(delta):
    start = 0.0
    end = 10.0
    filtrations = [start+i*delta for i in range(int((end-start)/delta) + 1)]
    return filtrations

def get_ref_eigs():
        filtrations = get_filtrations(delta=0.1)
        dims = [0, 1, 2]
        
        ref_eigs = [[],#dim 0
                        [],#dim 1
                        []]#dim 2
        for filtration in filtrations:
                for dim in dims:
                        filename = f"./protein_L/dim{dim}_a{filtration:.2f}_b{filtration+0.1:.2f}.mkt"

                        L = mmread(filename)
                        eigs = scipy.linalg.eigvalsh(L.todense())
                        ref_eigs[dim].append(eigs)
                        print(f"Got d={dim}, a={filtration}")
        return ref_eigs




def get_time_max_error_failures(tol, all_ref_eigs, min_ncv=200):
        filtrations = get_filtrations(delta=0.1)
        total_time = 0
        max_error = 0
        dims = [1,2]
        failures = 0
        for idx, filtration in enumerate(filtrations):
                for dim in dims:
                        # read file
                        filename = f"./protein_L/dim{dim}_a{filtration:.2f}_b{filtration+0.1:.2f}.mkt"
                        # print(f"Try to read file: {filename}")
                        L = mmread(filename)
                        # print(f"File {filename} has shape {L.shape}")
                        simplices = L.shape[0]
                        if simplices == 0 or L.count_nonzero() == 0:
                                continue
                        
                        k=10
                        k = min(k, simplices - 1) if simplices > 1 else 1
                        # ncv = 209
                        ncv = min(simplices, min_ncv)
                        print(f"do d={dim}, a={filtration}, k={k},tol={tol},ncv={ncv},simplices={simplices}")
                        t = time.time()
                        sigma = 1e-3
                        try:
                                eigs = scipy.sparse.linalg.eigsh(L, k=k, which="LM", ncv=ncv, sigma=sigma, tol=tol, maxiter=2000,return_eigenvectors=False)
                                elapsed_time = time.time() - t
                                print(f"Time taken: {elapsed_time}")
                                total_time = total_time + elapsed_time
                                
                                ref_eigs = all_ref_eigs[dim][idx][0:k]
                                residual = [abs(eigs[i] - ref_eigs[i]) for i in range(len(eigs))] 
                                
                                current_error = max(residual)
                                if current_error > 1e-1:
                                        print(f"Eigs: {eigs}")
                                        print(f"ref_eigs: {ref_eigs}")
                                        print(f"residual: {residual}")
                                if current_error > max_error:
                                        max_error = current_error
                                        print(f"new max_error: {max_error}")
                        except Exception as e:
                                elapsed_time = time.time() - t
                                print(f"Time taken: {elapsed_time}")
                                total_time = total_time + elapsed_time
                                failures = failures+1
                                print(e)    
                                # raise e            
                        
        return total_time, max_error, failures

def do_test(tol):
    all_ref_eigs = get_ref_eigs()
    data = []

    min_ncvs = [20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250]#, 300, 500]
    for min_ncv in min_ncvs:
        t, e, f = get_time_max_error_failures(tol, all_ref_eigs, min_ncv)
        data.append([tol, min_ncv, t, e,f ])
        print(f"Tol={tol}, t={t}, e={e}, min_ncv={min_ncv},failures={f}",flush=True)

    df = pd.DataFrame(data,columns=[str(tol),"min_ncv", "time", "max_error", "failures"])
    df.to_csv(f"profiles/tol{tol}.csv",index=False)


if __name__ == "__main__":
    tol = float(sys.argv[1]) # python ncv_tol_test.py <tol>
    do_test(tol=tol)    