import gudhi
import petls
from math import sqrt
import numpy as np
import pytest

def my_restriction(simplex: list[int], coface: list[int], sst: petls.sheaf_simplex_tree) -> float:
    from math import sqrt
    if len(simplex) == 1:
        if simplex == [coface[0]]:
            sibling = [coface[1]]
        else:
            sibling = [coface[0]]
        
        coords_simplex = sst.extra_data[tuple(simplex)][0:3]
        coords_sibling = sst.extra_data[tuple(sibling)][0:3]
        distance = sqrt((coords_simplex[0] - coords_sibling[0])**2 \
                    + (coords_simplex[1] - coords_sibling[2])**2 \
                    + (coords_simplex[2] - coords_sibling[1])**2)
        return sst.extra_data[tuple(sibling)][3] / distance
    elif len(simplex) == 2:
        coeff = 1.0
        for (sibling, _) in sst.st.get_boundaries(coface):
            if sibling == simplex:
                opposite_vertex = coface[sst.coface_index(simplex,coface)]
                coeff = coeff * sst.extra_data[tuple([opposite_vertex])][3] #charge
            else:
                coeff = coeff / sst.st.filtration(sibling)
        return coeff

    return 1

def get_sst(points, charges):
    st = gudhi.RipsComplex(points=points, max_edge_length=6).create_simplex_tree(max_dimension=3)


    extra_data = {
        tuple([0]): [*points[0],charges[0]],
        tuple([1]): [*points[1],charges[1]],
        tuple([2]): [*points[2],charges[2]],
        tuple([0,1]): 1,
        tuple([0,2]): 1,
        tuple([1,2]): 1,
        tuple([0,1,2]): 0
    }

    sst = petls.sheaf_simplex_tree(st,extra_data,my_restriction)
    return sst

def test_sst():
    points = [[0,0,0], [3, 0, 0], [0, 4, 0]]
    as_np = [np.array(x) for x in points]
    dists = [np.linalg.norm(as_np[0]-as_np[1]), np.linalg.norm(as_np[0]-as_np[2]), np.linalg.norm(as_np[1]-as_np[2])] 
    charges = [2, 7, 11]
    expected_cbdys = [
        np.array([[-charges[1]/dists[0], charges[0]/dists[0], 0],
        [-charges[2]/dists[1], 0, charges[0]/dists[1]],
        [0, -charges[2]/dists[2], charges[1]/dists[2]]
        ]),
        np.array([[charges[2]/(dists[1]*dists[2]), -charges[1]/(dists[0]*dists[2]), charges[0]/(dists[0]*dists[1])]])
    ]
    expected_filtrations = [[0,0,0], [3,4,5], [5]]
    sst = get_sst(points, charges)
    coboundaries, filtrations = sst.apply_restriction_function()
        
    for i in range(len(coboundaries)):
        np.testing.assert_allclose(coboundaries[i], expected_cbdys[i])
    for i in range(len(filtrations)):   
        np.testing.assert_allclose(np.array(filtrations[i]), np.array(expected_filtrations[i]))

if __name__ == "__main__":
    test_sst()