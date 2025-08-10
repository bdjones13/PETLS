from petls import Rips
import numpy as np
import scipy
import pytest

def compare_spectra(ref, test):
    print(f"compare {ref} to {test}") # only actually prints if assert fails in pytest
    if ref[0] != test[0] or ref[1] != test[1] or ref[2] != test[2]:
        print("failed.")
        return False
    assert ref[3] == pytest.approx(test[3])

def compare_spectra_multiple(ref, test):
    for i in range(len(ref)):
        compare_spectra(ref[i],test[i])


def rips_points(threshold = None):
    max_dim = 3
    points = np.array([
        [0, 0],
        [0, 3],
        [4, 0],
        [4, 3]])
    return Rips(points=points, max_dim=max_dim, threshold=threshold)

def rips_distances(threshold = None):
    max_dim = 3
    distances = np.array(
        [[0,0,0,0],
        [3,0,0,0],
        [4,5,0,0],
        [5,4,3,0]]
    )
    return Rips(distances = distances, max_dim=max_dim, threshold=threshold)

def rips_file(threshold = None):
    import os
    current_directory = os.getcwd()
    # print(current_directory)
    max_dim = 3
    return Rips(filename="data/rips/rect.lower_distance_matrix", max_dim=max_dim, threshold=threshold)

def pl_compare(pl):
    expected = [
        (0, 0, 3, [0,0,2,2]),
        (1, 0, 3, []),
        (2, 0, 3, []),
        (3, 0, 3, []),
        (0, 3, 4, [0,2,2,4]),
        (1, 3, 4, [2,2]),
        (2, 3, 4, []),
        (3, 3, 4, []),
        (0, 4, 5, [0,4,4,4]),
        (1, 4, 5, [2,2,4,4]),
        (2, 4, 5, []),
        (3, 4, 5, []),
        (0, 5, 5, [0,4,4,4]),
        (1, 5, 5, [4,4,4,4,4,4]),
        (2, 5, 5, [4,4,4,4]),
        (3, 5, 5, [4]),
    ]
    compare_spectra_multiple(expected, pl.spectra())

def pl_compare_threshold(pl):
    expected = [
        (0, 0, 3, [0,0,2,2]),
        (1, 0, 3, []),
        (2, 0, 3, []),
        (3, 0, 3, []),
        (0, 3, 4, [0,2,2,4]),
        (1, 3, 4, [2,2]),
        (2, 3, 4, []),
        (3, 3, 4, []),
        (0, 4, 4, [0,2,2,4]),
        (1, 4, 4, [0,2,2,4]),
        (2, 4, 4, []),
        (3, 4, 4, []),
    ]
    compare_spectra_multiple(expected, pl.spectra())

def test_rips():
    pls = [rips_points(), rips_distances(), rips_file()]
    for pl in pls:
        pl_compare(pl)
    
    threshold = 4.5
    pls = [rips_points(threshold), rips_distances(threshold), rips_file(threshold)]
    for pl in pls:
        pl_compare_threshold(pl)

    

if __name__ == "__main__":
    test_rips()