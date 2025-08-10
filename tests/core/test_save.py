import petls
import numpy as np
import scipy

d1 = np.array([[-1,0,-1],
               [1,-1,0],
               [0,1,1]])
d2 = np.array([[1],[1],[-1]]) 
boundaries = [d1,d2]
filtrations = [[0,1,2],[3,4,5],[5]]

pl = petls.Complex(boundaries, filtrations)
pl.store_L(0,0,1,"saved_matrix.mtx")