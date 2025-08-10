import petls
import numpy as np
import scipy
import pytest
import os
from pathlib import Path

def get_complex():
    d1 = np.array([[-1,0,-1],
                [1,-1,0],
                [0,1,1]])
    d2 = np.array([[1],[1],[-1]]) 
    boundaries = [d1,d2]
    filtrations = [[0,1,2],[3,4,5],[5]]
    complex = petls.Complex(boundaries, filtrations)
    return complex

def test_profile():
    complex = get_complex()
    complex.spectra(0,1,2)
    complex.spectra(0,2,3)
    complex.spectra(1,2,3)

    filename = "test_profile.csv"
    if os.path.exists(filename):
        os.remove(filename)

    complex.time_to_csv(filename)

    assert os.path.exists(filename)
    # if it exists, remove it
    os.remove(filename)

if __name__ == "__main__":
    test_profile()