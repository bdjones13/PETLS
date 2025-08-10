import petls
import numpy as np
import scipy
import pytest

def compare_spectra(ref, test):
    print(f"compare {ref} to {test}") # only actually prints if assert fails in pytest or if run manually
    if ref[0] != test[0] or ref[1] != test[1] or ref[2] != test[2]:
        return False
    assert ref[3] == pytest.approx(test[3])

def compare_spectra_multiple(ref, test):
    for i in range(len(ref)):
        compare_spectra(ref[i],test[i])


def get_pl():
    d1 = np.array([[-1,0,-1],
                [1,-1,0],
                [0,1,1]])
    d2 = np.array([[1],[1],[-1]]) 
    boundaries = [d1,d2]
    filtrations = [[0,1,2],[3,4,5],[5]]
    pl = petls.Complex(boundaries, filtrations)
    return pl

def test_specific():
    pl = get_pl()
    assert pl.spectra(0,5,6) == pytest.approx([0,3,3])

def test_all():
    # print("test_all")
    pl = get_pl()
    # print("got pl")
    
    # print("print pl boundaries:",flush=True)
    # pl.print_boundaries()
    # print("finished printing boundaries???",flush=True)
    s = pl.spectra()
    # print("got spectra")
    # print("spectra: ", s)
    expected = [
        (0, 0, 3, [0]),
        (1, 0, 3, []),
        (2, 0, 3, []),
        (0, 3, 4, [0, 1, 3]),
        (1, 3, 4, [2]),
        (2, 3, 4, []),
        (0, 4, 5, [0, 3, 3]),
        (1, 4, 5, [1, 3]),
        (2, 4, 5, []),
        (0, 5, 5, [0, 3, 3]),
        (1, 5, 5, [3, 3, 3]),
        (2, 5, 5, [3])
    ]
    compare_spectra_multiple(expected, s)

def test_L():
    pl = get_pl()
    L = pl.get_L(1,5,5)
    ref = np.diag(3*np.ones(3)) # diagonal 3x3 with entries all 3
    np.testing.assert_allclose(L,ref)

def test_up():
    pl = get_pl()
    up = pl.get_up(0, 1, 3)
    # print("up(1,5,5)=",pl.get_up(1,5,5))
    ref = np.array([[1,-1], [-1, 1]])
    np.testing.assert_allclose(up,ref)

def test_down():
    pl = get_pl()
    down = pl.get_down(1,5)
    # print("down=",down)
    ref = np.array([[2, -1, 1],
                       [-1, 2, 1],
                       [1, 1, 2]])
    np.testing.assert_allclose(down, ref)

def test_sum_up_down():
    pl = get_pl()
    down = pl.get_down(1,5)
    up = pl.get_up(1,5,5)
    L = pl.get_L(1,5,5)
    np.testing.assert_allclose(L, up+down)

def test_eigenvectors():
    pl = get_pl()
    x = pl.eigenpairs()
    for pair in x:
        dim = pair[0]
        a = pair[1]
        b = pair[2]
        eigs = pair[3]
        eigenvectors = pair[4]
        # print("a: ", type(a),a)
        # print("b: ", type(b),b)
        # print("dim: ", type(dim),dim)
        # print("eigs: ", type(eigs),eigs)
        # print("eigenvectors: ", type(eigenvectors),eigenvectors)

def test_allpairs():
    pl = get_pl()
    s = pl.spectra(allpairs=True)
    expected = [(0, 0, 0, [0]), 
                (1, 0, 0, []), 
                (2, 0, 0, []), 
                (0, 0, 3, [0]), 
                (1, 0, 3, []), 
                (2, 0, 3, []),
                (0, 0, 4, [0]),
                (1, 0, 4, []),
                (2, 0, 4, []),
                (0, 0, 5, [0]),
                (1, 0, 5, []),
                (2, 0, 5, []), 
                (0, 3, 3, [0, 0, 2]),
                (1, 3, 3, [2]),
                (2, 3, 3, []),
                (0, 3, 4, [0, 1, 3]),
                (1, 3, 4, [2]), 
                (2, 3, 4, []), 
                (0, 3, 5, [0, 3, 3]), 
                (1, 3, 5, [2]), 
                (2, 3, 5, []),
                (0, 4, 4, [0, 1, 3]),
                (1, 4, 4, [1, 3]),
                (2, 4, 4, []), 
                (0, 4, 5, [0, 3, 3]),
                (1, 4, 5, [1, 3]),
                (2, 4, 5, []), 
                (0, 5, 5, [0, 3, 3]), 
                (1, 5, 5, [3, 3, 3]), 
                (2, 5, 5, [3])]
    compare_spectra_multiple(expected, s)

if __name__ == "__main__":
    # test_allpairs()
    # test_all()
    # test_L()
    # test_up()
    test_down()
    # test_sum_up_down()
    # test_eigenvectors()