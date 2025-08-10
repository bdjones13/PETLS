from petls import sheaf_simplex_tree, NewPersistentSheafLaplacian
from test_sst import get_sst
import pytest
import numpy as np
import gudhi
from math import sqrt

def test_psl():
    points = [[0,0,0], [3, 0, 0], [0, 4, 0]]
    charges = [2, 7, 11]
    
    as_np = [np.array(x) for x in points]
    dists = [np.linalg.norm(as_np[0]-as_np[1]), np.linalg.norm(as_np[0]-as_np[2]), np.linalg.norm(as_np[1]-as_np[2])] 
    

    sst = get_sst(points, charges) # see test_sst.py for details    

    psl = NewPersistentSheafLaplacian(sst)
    
    q0 = charges[0]
    q1 = charges[1]
    q2 = charges[2]
    d01 = dists[0]
    d02 = dists[1]
    d12 = dists[2]
    
    cbdy0 = np.array([[-q1/d01, q0/d01, 0],
    [-q2/d02, 0, q0/d02],
    [0, -q2/d12, q1/d12]])
    cbdy1 = np.array([[q2/(d02*d12), -q1/(d01*d12), q0/(d01*d02)]])
    
    bdy1 = cbdy0.T
    bdy2 = cbdy1.T

    expected_L0 = np.matmul(bdy1,bdy1.T)
    expected_L1 = bdy1.T @ bdy1 + bdy2 @ bdy2.T
    expected_L2 = bdy2.T @ bdy2    
    
    np.testing.assert_allclose(psl.get_L(0,5,5), expected_L0)
    np.testing.assert_allclose(psl.get_L(1,5,5), expected_L1)
    np.testing.assert_allclose(psl.get_L(2,5,5), expected_L2)


if __name__ == "__main__":
    test_psl()