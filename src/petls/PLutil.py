import gudhi.simplex_tree
from scipy.sparse import coo_matrix
import pytest
import numpy as np
import matplotlib.pyplot as plt

def summaries(spectra, func):
    """ Apply a function to the eigenvalues of L_{dim}^{a,b} for all (dim, a, b, eigs) in spectra.

        Parameters
        ----------
        spectra : List[Tuple[int, float, float, List[float]]]
            List of (dim, a, b, eigs) to apply func.
        func : Callable
            Function to apply to the eigenvalues of L_{dim}^{a,b} o for every tuple (dim, a, b) in spectra.

        Returns
        -------
        List[np.array]
            List of numpy arrays, where row b, column a is the result of applying func to the eigenvalues of L_{dim}^{a,b}.
        List[float]
            List of all filtrations in the input spectra. This is an ordered index of the summary arrays.
        List[int]
            List of dimensions in the input spectra. There is one summary array for each dimension.

        Example
        -------
        >>> def min_nonzeros(eigs):
                tol = 1e-6
                nonzeros = [x for x in eigs if x > tol]
                if len(nonzeros) == 0:
                    return 0
                return min(nonzeros)
        >>> spectra = [(0, 0, 0, [0]), (0, 0, 3, [0]), (0, 0, 4, [0]), (0, 3, 3, [0, 0, 2]), (0, 3, 4, [0, 1, 3]), (0, 4, 4, [0, 1, 3])]
        >>> summaries(spectra, min_nonzeros)
        ([array([[0.        , 0.        , 0.        ],
                [0.        , 1.99999988, 0.        ],
                [0.        , 0.99999988, 0.99999988]])], [0.0, 3.0, 4.0], [0])
        """
    all_filtrations = set()
    dims = set()
    for dim, a, b, eigs in spectra:
        all_filtrations.add(a)
        all_filtrations.add(b)
        dims.add(dim)
    all_filtrations = sorted(list(all_filtrations))
    indexed_filtrations = {a: i for i, a in enumerate(all_filtrations)} # dict with filtration as key and index as value
    indexed_dims = {dim: i for i, dim in enumerate(sorted(list(dims)))}
    num_filtrations = len(all_filtrations)
    summaries = [np.zeros((num_filtrations,num_filtrations)) for _ in range(len(dims))]

    for dim, a, b, eigs in spectra:
        summaries[indexed_dims[dim]][indexed_filtrations[b], indexed_filtrations[a]] += func(eigs) # addition modifies actual matrix, not just a view of it
    
    return summaries, all_filtrations, sorted(dims)

def plot_summary(ax, summary, **kwargs):
    """ Wrapper for matplotlib imshow to plot a summary of the spectra in the shape of a persistence diagram with diagonal displayed.
        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes to plot the summary on.
        summary : np.array
            Summary of the spectra in the shape of a persistence diagram. Row b, column a is the summary of L^{a,b}.
        **kwargs : keyword arguments
            Additional keyword arguments to pass to imshow. For example, cmap="Blues" to set the colormap.
        Returns
        ------- 
        pos : matplotlib.image.AxesImage
            The imshow object. This can be used to set the colorbar or other properties of the plot.

    """
    pos = ax.imshow(summary, origin="lower", **kwargs)
    ax.plot([0,1],[0,1], transform=ax.transAxes,color="black") # plot y=x
    return pos


def simplex_tree_boundaries_filtrations(simplex_tree: gudhi.simplex_tree):
    """ Helper function for extracting the boundary matrix from a Gudhi simplex tree. 
    This is used only by the Complex constructor, and is not a function the user needs access to.
    
    """
    # give each simplex a unique index
    index = 0
    indices = {}
    for simplex_with_filtration in simplex_tree.get_filtration():
        indices[tuple(simplex_with_filtration[0])] = index
        index = index + 1
    # set up with space for boundaries and filtrations in each dimension 
    boundaries_triples = [[] for _ in range(simplex_tree.dimension()+1)]
    filtrations = [[] for _ in range(simplex_tree.dimension()+1)]    

    # get all boundaries in format (row, col, coeff)
    for simplex_with_filtration in simplex_tree.get_filtration():
        simplex = simplex_with_filtration[0]
        filtration = simplex_with_filtration[1]
        dim = len(simplex) - 1
        filtrations[dim].append(filtration)

        if dim == 0: # no boundaries from d_0
            continue

        coeff = 1
        for face, filtration in simplex_tree.get_boundaries(simplex):
            # print(dim, indices[tuple(face)], indices[tuple(simplex)])
            # print(len(boundaries_triples))
            boundaries_triples[dim].append([indices[tuple(face)], indices[tuple(simplex)], coeff])            
            coeff *= -1

    # now create dictionaries that map actually-occuring simplex indices bijectively to [0, 1, ..., N]  
    index_mappings = [{} for _ in range(simplex_tree.dimension() + 1)]

    for dim in range(simplex_tree.dimension()):
        indices_of_actual_simplices = list(set([triple[0] for triple in boundaries_triples[dim+1]]).union(
                                           set([triple[1] for triple in boundaries_triples[dim]])))
        for i in range(len(indices_of_actual_simplices)):
            index_mappings[dim][indices_of_actual_simplices[i]] = i
    
    indices_of_actual_simplices = list(set([triple[1] for triple in boundaries_triples[simplex_tree.dimension()]]))
    for i in range(len(indices_of_actual_simplices)):
        index_mappings[simplex_tree.dimension()][indices_of_actual_simplices[i]] = i

    # now convert to coo format (lists for row, col, data) then numpy array
    boundaries = []
    for dim in range(1,simplex_tree.dimension()+1):
        row = []
        col = []
        data = []
        boundary_triples = boundaries_triples[dim]
        for triple in boundary_triples:
            row.append(index_mappings[dim-1][triple[0]])
            col.append(index_mappings[dim][triple[1]])
            data.append(triple[2])
        boundary = coo_matrix((data, (row,col)), shape=(len(index_mappings[dim-1]),
                                                              len(index_mappings[dim])
                                                                )).toarray()
        boundaries.append(boundary)
    return boundaries, filtrations


def compare_spectra(ref, test):
    print(f"compare {ref} to {test}") # only actually prints if assert fails in pytest
    if ref[0] != test[0] or ref[1] != test[1] or ref[2] != test[2]:
        print("failed.")
        return False
    assert ref[3] == pytest.approx(test[3])

def compare_spectra_multiple(ref, test):
    for i in range(len(ref)):
        compare_spectra(ref[i],test[i])