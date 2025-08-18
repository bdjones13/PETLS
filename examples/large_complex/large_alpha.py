
import petls
import tadasets
import scipy
import numpy as np


def get_complex(n):
    replicate = 0
    complex = petls.Alpha(torus, max_dim=3)
    return complex

def get_filtrations(delta):
    start = 0.0
    end = 4.0 #5.0
    filtrations = [start+i*delta for i in range(int((end-start)/delta) + 1)]
    return filtrations

def store_eigs(all_eigs, filename):
    with open(filename, 'w') as f:
        for s in all_eigs:
            f.write(f"{s[0]},{s[1]},{s[2]},{','.join(map(str, s[3]))}\n")



n = 500

torus = tadasets.torus(n=n, c=3, a=1, noise=0.0, seed=0)
complex = petls.Alpha(torus, max_dim=3)

filtrations = get_filtrations(delta=0.25)
complex.set_eigs_Algorithm(lambda L: scipy.linalg.eigvalsh(L))
dims = [0, 1, 2]
# compute all the eigenvalues
all_eigs = []
for idx, filtration in enumerate(filtrations):
    for idx2, filtration2 in enumerate(filtrations):
        if (idx2 < idx) or (idx == 0 and idx2 == 0):
            continue
        for dim in dims:      
            s = complex.spectra(dim, a=filtration, b=filtration2)
            all_eigs.append(s)

store_eigs(all_eigs, 'all_eigs.txt')