from petls import Complex
import numpy as np

d1 = np.array([[-1,  0, -1],
                [1, -1,  0],
                [0,  1,  1]])

d2 = np.array([[1], [1], [-1]]) 
boundaries = [d1,d2]
filtrations = [[0,1,2], # c0
               [3,4,5], # c1
               [5]]     # c2

pl = Complex(boundaries, filtrations)
eigenvalues = pl.spectra(dim=1, a=4, b=5)
print(eigenvalues)