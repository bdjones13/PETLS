from petls import Alpha
import numpy as np

d1 = np.array([[-1,0,-1],
               [1,-1,0],
               [0,1,1]])
d2 = np.array([[1],[1],[-1]]) 
boundaries = [d1,d2]
filtrations = [[0,1,2],[3,4,5],[5]]

complex = Alpha(filename="data/alpha/input",max_dim=3)
print(complex.spectra())

points = [[0.0, 0.0, 0.0],
[0.0, 0.0, 1.0],
[0.0, 1.0, 0.0],
[1.0, 0.0, 0.0]]

complex_points = Alpha(points=points,max_dim=3)
complex_points.verbose = True
print(complex_points.spectra())