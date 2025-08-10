from ._petls import *

import enum

########################
#  for sst
import gudhi.simplex_tree
from scipy.sparse import coo_matrix
from collections.abc import Callable # for the type hint of the restriction function
import numpy as np
########################
from .PLutil import simplex_tree_boundaries_filtrations, summaries, plot_summary

# Define the behavior of "import petls" and options for "from petls import X"
__all__ = ["Complex", "Rips", "Alpha", "dFlag", "sheaf_simplex_tree", "PersistentSheafLaplacian","NewPersistentSheafLaplacian", "summaries", "plot_summary"]


# The petls library is written primarily in C++.
# A layer of Python bindings (via pybind11) is used to create many Python classes that correspond
# to multiple C++ template options.
# Then a second layer of Python is added to make the first layer user-friendly. This second layer is what this file contains.   
# The "first" Python layer is defined in  src/_petls.cpp, src/core/Complex.cpp, and src/variants/*.


class eigs_Algorithms(enum.Enum):
    """
    Enum to choose which eigenvalue algorithm to use.
    """
    selfadjoint = 1
    eigensolver = 2
    bdcsvd = 3

class up_Algorithms(enum.Enum):
    """
    Enum to choose which up-Laplacian algorithm to use
    """
    schur = 1


class storage(enum.Enum):
    int = 1
    float = 2


class Complex(object):
    """Primary class used to compute persistent Laplacian matrices and eigenvalues from a complex.
    
    Attributes
    ----------
    verbose : boolean
        Print progress if spectra() is called
    flipped : boolean
        Compute the top-dimensional Laplacian's eigenvalues via the eigenvalues of the smaller of B_N B_N^T or B_N^T B_N and possible zero-padding
        
    Methods
    -------
    set_boundaries_filtrations(boundaries, filtrations)
        If the boundaries and filtrations were not set in the constructor, set them here.
    get_L(dim, a, b)
        Get the persistent Laplacian matrix.
    get_up(dim, a, b)
        Get the persistent up-Laplacian matrix.
    get_down(dim, a)
        Get the persistent down-Laplacian matrix.
    nonzero_spectra(dim, a, b, PH_basis=None, use_dummy_harmonic_basis=True)
        Compute the nonzero eigenvalues of the PL using reduction by persistent homology or standard basis of the null space.
    spectra(dim = None, a = None, b = None, request_list = None)
        Compute the eigenvalues of L_{dim}^{a,b} or for every tuple (dim, a, b) in request_list.
    eigenpairs(dim = None, a = None, b = None, request_list = None)
        Compute the eigenvalues and eigenvectors of L_{dim}^{a,b} or for every tuple (dim, a, b) in request_list.
    eigenvalues_summarize(eigenvalues)
        Compute the betti number and least nonzero eigenvalue of a list of eigenvalues.
    print_boundaries(self):
        Print all boundaries and corresponding filtrations in the complex.
    store_spectra(spectra_list, file_prefix):
        Store all of the eigenvalues in files f"{file_prefix}_spectra_{dim}.txt" for each dimension in the complex.    
    store_spectra_summary(spectra_list, file_prefix):
        Store eigenvalue summaries in file f"{file_prefix}_spectra_summary.txt".
    time_to_csv(self, filename):
        Store the time taken to compute the spectra in a CSV file. This is done automatically when spectra() is called.
    filtration_list_to_spectra_request(filtrations, dims):
        Get a list of tuples (dim, a, b) for all combinations of dimension and successive filtration values (a=filtrations[i], b=filtrations[i+1]).
    get_all_filtrations(self):
        Get a sorted list of all filtration values that occur in the complex.
    """


    def __init__(self, boundaries = None, filtrations = None, eigs_Algorithm = eigs_Algorithms.selfadjoint, up_Algorithm = up_Algorithms.schur, simplex_tree = None, storage = None):
        if storage is not None:
            pl_class = getattr(_petls, "Complex")
        else:
            pl_class = getattr(_petls, "Complex")
        
        if boundaries is None and filtrations is None:
            self.pl = pl_class()
        else:
            self.pl = pl_class(boundaries, filtrations)
        self.verbose = False
        self.flipped = False

        if simplex_tree is not None:
            boundaries, filtrations = simplex_tree_boundaries_filtrations(simplex_tree)
            self.set_boundaries_filtrations(boundaries, filtrations)

    def set_boundaries_filtrations(self, boundaries, filtrations):
        """ If the boundaries and filtrations were not set in the constructor, set them here.

        Parameters
        ----------
        boundaries : List[np.array]
            List of boundary matrices
        filtrations : List[List[float]]
            For each dimension in the complex, a list of simplex filtration values

        Returns
        -------
        None

        Examples
        --------
        >>> d1 = np.array([[-1,0,-1],
                [1,-1,0],
                [0,1,1]])
        >>> d2 = np.array([[1],[1],[-1]]) 
        >>> boundaries = [d1,d2]
        >>> filtrations = [[0,1,2],[3,4,5],[5]]
        >>> pl = petls.Complex()
        >>> pl.set_boundaries_filtrations(boundaries, filtrations)

        """
        self.pl.set_boundaries_filtrations(boundaries, filtrations)
    
    @property
    def verbose(self):
        return self._verbose

    @verbose.setter
    def verbose(self, verbose):
        self._verbose = verbose
        self.pl.set_verbose(verbose)

    @property
    def flipped(self):
        return self._flipped
    
    @flipped.setter
    def flipped(self, flipped):
        self._flipped = flipped
        self.pl.set_flipped(flipped)

    def set_flipped(self, flipped=True):
        self.pl.set_flipped(flipped)

    def get_L(self, dim, a, b):
        return self.pl.get_L(dim, a, b)
    
    def get_up(self, dim, a, b):
        return self.pl.get_up(dim, a, b)

    def get_down(self, dim, a):
        return self.pl.get_down(dim, a)

    def nonzero_spectra(self, dim, a, b, PH_basis=None, use_dummy_harmonic_basis=True):
        self.pl.nonzero_spectra(dim, a, b, PH_basis, use_dummy_harmonic_basis)


    def spectra(self,dim = None, a = None, b = None, request_list = None, allpairs = False):
        """ Compute the eigenvalues of L_{dim}^{a,b} or for every tuple (dim, a, b) in request_list.

        Parameters
        ----------
        dim : int, optional
            Dimension.
        a : float, optional
            Start filtration value.
        b : float, optional
            End filtration value
        request_list : List[List[int, float, float]], optional
            List of (dim, a, b) to compute the spectra of

        Returns
        -------
        List[float]
            If passed dim, a, and b, returns eigenvalues of L_{dim}^{a,b}
        List[Tuple[int, float, float, List[float]]]
            If passed request_list, returns a tuple (dim, a, b, eigenvalues) for each request in request_list    
            If passed no arguments, returns as if request_list is of all combinations (dim, a, b) where b is the next filtration after a, and (dim, a, a) when a is the largest filtration value  
        
        Examples
        --------
        >>> d1 = np.array([[-1,0,-1],
                [1,-1,0],
                [0,1,1]])
        >>> d2 = np.array([[1],[1],[-1]]) 
        >>> boundaries = [d1,d2]
        >>> filtrations = [[0,1,2],[3,4,5],[5]]
        >>> pl = Complex(boundaries, filtrations)
        >>> pl.spectra(0, 1.2, 4.5)
        [0.0,  1.9999998807907104]
        >>> pl.spectra(request_list = [[0, 1.2, 4.5]])
        [(0, 1.2, 4.5, [0.0, 1.9999998807907104])]
        >>> pl.spectra()
        [(0, 0.0, 3.0, [0.0]), (1, 0.0, 3.0, []), (2, 0.0, 3.0, []), (0, 3.0, 4.0, [0.0, 0.9999998807907104, 3.0]), (1, 3.0, 4.0, [2.0]), (2, 3.0, 4.0, []), (0, 4.0, 5.0, [0.0, 2.999999761581421, 3.0]), (1, 4.0, 5.0, [0.9999999403953552, 2.999999761581421]), (2, 4.0, 5.0, []), (0, 5.0, 5.0, [0.0, 2.999999761581421, 3.0]), (1, 5.0, 5.0, [3.0, 3.0, 3.0]), (2, 5.0, 5.0, [3.0])]
        """
        # print("__init__.py pl.spectra",flush=True)
        
        if request_list is not None:
            # print("request_list is not none")
            return self.pl.spectra(request_list)
        elif dim is not None and a is not None and b is not None:
            # print("dim, a, b is not none")
            return self.pl.spectra(dim, a, b)
        elif allpairs:
            return self.pl.spectra_allpairs()
        else:
            # print("call self.pl.spectra()")
            return self.pl.spectra()
        
    def eigenpairs(self,dim = None, a = None, b = None, request_list = None, allpairs = False):
        """ Compute the eigenvalues and eigenvectors of L_{dim}^{a,b} or for every tuple (dim, a, b) in request_list.

        Parameters
        ----------
        dim : int, optional
            Dimension.
        a : float, optional
            Start filtration value.
        b : float, optional
            End filtration value
        request_list : List[List[int, float, float]], optional
            List of (dim, a, b) to compute the spectra of

        Returns
        -------
        Tuple[List[float], numpy.ndarray]
            If passed dim, a, and b, returns eigenvalues and eigenvectors of L_{dim}^{a,b}
        List[Tuple[int, float, float, List[float], numpy.ndarray]]
            If passed request_list, returns a tuple (dim, a, b, eigenvalues, eigenvectors) for each request in request_list.   
            If passed no arguments, returns as if request_list is of all combinations (dim, a, b) where b is the next filtration after a, and (dim, a, a) when a is the largest filtration value.  
            If passed only one argument, assume it is a request_list.
        
        Examples
        --------
        >>> d1 = np.array([[-1,0,-1],
                [1,-1,0],
                [0,1,1]])
        >>> d2 = np.array([[1],[1],[-1]]) 
        >>> boundaries = [d1,d2]
        >>> filtrations = [[0,1,2],[3,4,5],[5]]
        >>> pl = Complex(boundaries, filtrations)
        >>> pl.spectra(0, 1.2, 4.5)
        [0.0,  1.9999998807907104]
        >>> pl.spectra(request_list = [[0, 1.2, 4.5]])
        [(0, 1.2, 4.5, [0.0, 1.9999998807907104])]
        >>> pl.spectra()
        [(0, 0.0, 3.0, [0.0]), (1, 0.0, 3.0, []), (2, 0.0, 3.0, []), (0, 3.0, 4.0, [0.0, 0.9999998807907104, 3.0]), (1, 3.0, 4.0, [2.0]), (2, 3.0, 4.0, []), (0, 4.0, 5.0, [0.0, 2.999999761581421, 3.0]), (1, 4.0, 5.0, [0.9999999403953552, 2.999999761581421]), (2, 4.0, 5.0, []), (0, 5.0, 5.0, [0.0, 2.999999761581421, 3.0]), (1, 5.0, 5.0, [3.0, 3.0, 3.0]), (2, 5.0, 5.0, [3.0])]
        """
        
        if (dim is not None and a is None and b is None and request_list is None):
            request_list = dim
            dim = None
        
        if request_list is not None:
            return self.pl.eigenpairs(request_list)
        elif dim is not None and a is not None and b is not None:
            return self.pl.eigenpairs(dim, a, b)
        else:
            return self.pl.eigenpairs()
        
        
    def eigenvalues_summarize(self, eigenvalues):
        return self.pl.eigenvalues_summarize(eigenvalues)
    
    def store_L(self, dim, a, b, prefix):
        # caution: recomputes L
        self.pl.store_L(dim, a, b, prefix)
    
    def print_boundaries(self):
        # print("in py print_boundaries",flush=True)
        self.pl.print_boundaries()
        # print("end py print_boundaries",flush=True)
    
    def store_spectra(self, spectra_list, file_prefix):
        self.pl.store_spectra(spectra_list, file_prefix)
    
    def store_spectra_summary(self, spectra_list, file_prefix):
        self.pl.store_spectra_summary(spectra_list, file_prefix)

    def time_to_csv(self, filename):
        self.pl.time_to_csv(filename)

    def filtration_list_to_spectra_request(self, filtrations, dims, allpairs = False):
        if allpairs:
            return self.pl.filtration_list_to_spectra_request_allpairs(filtrations, dims)
        return self.pl.filtration_list_to_spectra_request(filtrations, dims)
    
    def get_all_filtrations(self):
        return self.pl.get_all_filtrations()
    


class dFlag(Complex):
    def __init__(self, filename, max_dim, eigs_Algorithm = eigs_Algorithms.selfadjoint, up_Algorithm = up_Algorithms.schur):
        pl_class = getattr(_petls, "dFlag")
        self.pl = pl_class(filename, max_dim)

class Rips(Complex):
    def __init__(self, filename = None, points = None, distances = None, max_dim=3, threshold = None, eigs_Algorithm = eigs_Algorithms.selfadjoint, up_Algorithm = up_Algorithms.schur):
        pl_class = getattr(_petls, "Rips")
        
        if threshold is None:
            if filename is not None:
                self.pl = pl_class(filename, max_dim)
            elif points is not None:
                self.pl = pl_class(points, max_dim)
            elif distances is not None:
                self.pl = pl_class.from_distances(distances, max_dim)
            else:
                raise ValueError('Rips requires either filename or point set as input')
        else:
            if filename is not None:
                self.pl = pl_class(filename, max_dim, threshold)
            elif points is not None:
                self.pl = pl_class(points,max_dim, threshold)            
            elif distances is not None:
                self.pl = pl_class.from_distances(distances, max_dim, threshold)
            else:
                raise ValueError('Rips requires either filename or point set as input')
            
class Alpha(Complex):
    def __init__(self, filename = None, points = None, max_dim=3):
        pl_class = getattr(_petls, "Alpha")
        if filename is not None:
            self.pl = pl_class(filename, max_dim)
        elif points is not None:
            self.pl = pl_class(points, max_dim)
        else:
            raise ValueError("Alpha complex requires filename or point set as input")


class sheaf_simplex_tree():
    """Wrap a Gudhi simplex tree with possibly added data and a given restriction function.

    Attributes
    ----------
    st : gudhi.simplex_tree
        Underlying simplex tree
    extra_data : dict
        Keys are simplices (list of int) converted to tuples, e.g. tuple([0,1,2]). Values can be anything.
    restriction : Callable[[list[int], list[int], 'sheaf_simplex_tree'], float]
        Restriction function from a simplex to a coface, which may need knowledge from the whole sheaf_simplex_tree. e.g. my_restriction(simplex, coface, sst)
    complex_dim : int
        Dimension of the complex, same as dimension of the simplex tree

    Methods
    -------
    coface_index(simplex, coface):
        Get the index of the missing vertex, e.g. coface_index([0,1,3],[0,1,2,3]) = 2.
    apply_restriction_function():
        Get coboundaries and filtrations
    """

    def __init__(self,_st: 'gudhi.simplex_tree', _extra_data: dict, _restriction: Callable[[list[int],list[int], 'sheaf_simplex_tree'], float]):
        self.st = _st
        self.extra_data = _extra_data
        self.restriction = _restriction
        self.complex_dim = _st.dimension()

        # give each simplex a unique index
        index = 0
        indices = {}
        for simplex_with_filtration in self.st.get_filtration():
            indices[tuple(simplex_with_filtration[0])] = index
            index = index + 1
        self.indices = indices
    
    def coface_index(self,simplex, coface):
        """Get the index of the missing vertex, e.g. coface_index([0,1,3],[0,1,2,3]) = 2.
        
        Parameters
        ----------
        simplex : list[int]
            simplex
        coface : list[int]
            coface
        
        Returns
        -------
        int
            Index of the missing vertex
        """

        if len(simplex) != len(coface)-1:
            raise ValueError(f"len(simplex) != len(coface)-1. len(simplex) = {len(simplex)}, len(coface)={len(coface)}")
        for i in range(len(simplex)):
            if simplex[i] != coface[i]:
                return i
        # last index
        return len(simplex)

    def apply_restriction_function(self) -> tuple[list[np.array],list[list[float]]] :
        """Get coboundaries and filtrations

        Returns
        -------
        tuple[list[np.array],list[list[float]]]
            The first element of the tuple is a list of coboundary matrices, the second element is a list of filtrations (one list[float] per dimension)

        Examples
        --------
        >>> st = gudhi.simplex_tree(...)
        >>> def my_restriction(simplex, coface, sst):
                ...
                return coeff # float
        >>> extra_data = {}
        >>> sst = sheaf_simplex_tree(st, extra_data, my_restriction)
        >>> [coboundaries, filtrations] = sst.apply_restriction_function()
        >>> print(coboundaries) 
            [   array([[-5., 0., 2.],
                       [0., -5., 3.],
                       [-3., 2., 0.]]), 
                array([[-3., 2., 5.]])
            ]
        >>> print(filtrations)
        [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [1.0]]

        """
        # get coboundaries

        # set up with space for coboundaries and filtrations in each dimension 
        coboundaries_triples = [[] for _ in range(self.complex_dim)]
        filtrations = [[] for _ in range(self.complex_dim+1)]

        # loop over all simplices
        for simplex_with_filtration in self.st.get_filtration():
            simplex = simplex_with_filtration[0]
            filtration = simplex_with_filtration[1]
            # print("simplex with filtration: (", simplex, ", ", filtration, ")")
            dim = len(simplex) - 1
            filtrations[dim].append(filtration)
            if dim == self.complex_dim: # no coboundary from top-dimension
                continue
            # loop over all cofaces
            for coface_with_filtration in self.st.get_cofaces(simplex,1):
                coface = coface_with_filtration[0]
                sign = (-1)**(self.coface_index(simplex, coface) % 2)

                # apply restriction function
                coeff = sign*self.restriction(simplex, coface, self)
                
                # store sparse matrix entry as [coface index, face index, value]
                coboundaries_triples[dim].append([self.indices[tuple(coface)], self.indices[tuple(simplex)], coeff])
        return self.reindex_coboundaries(coboundaries_triples), filtrations


    def reindex_coboundaries(self, coboundaries_triples):
        # This method is not documented with a formal docstring because it is meant to be a private method

        # Each simplex has a unique index with respect to the whole complex,
        # but we need to index the coboundary matrix with respect to each dimension's simplices
        # e.g. the matrix d2 may be [1, -1, 1] stored as [[5, 2, 1], [5, 3, -1], [5, 4, 1]],
        # this will convert it to [[0, 0, 1], [0, 1, -1 ], [0, 2, 1]],
        # then return as dense numpy array

        # first get all simplices in each dimension. It is possible that
        # an n-simplex does not have any (n+1)-simplex cofaces, so it would not appear
        # in any of the dimension n coboundary triples (from n to n+1), but it may be
        # the coface of an (n-1)-simplex and appear in the (n-1)-coboundary triples (from n-1 to n).
        # so we add all of them to a set.
        indices_of_actual_simplices_set = [set() for _ in range(self.complex_dim+1)]
        for dim in range(self.complex_dim):
            coboundary_triples = coboundaries_triples[dim]
            for triple in coboundary_triples:
                indices_of_actual_simplices_set[dim].add(triple[1])
                indices_of_actual_simplices_set[dim+1].add(triple[0])
        indices_of_actual_simplices = [list(indices) for indices in indices_of_actual_simplices_set]
        
        # now create dictionaries that map the simplex indices bijectively to [0, 1, ..., N]  
        index_mappings = [{} for _ in range(self.complex_dim + 1)]
        for dim in range(self.complex_dim+1):
            indices = indices_of_actual_simplices[dim]
            for i in range(len(indices)):
                index_mappings[dim][indices[i]] = i
        
        # convert to coo format (lists for row, col, data)
        coboundaries = []
        for dim in range(self.complex_dim):
            row = []
            col = []
            data = []
            coboundary_triples = coboundaries_triples[dim]
            for triple in coboundary_triples:
                row.append(index_mappings[dim+1][triple[0]])
                col.append(index_mappings[dim][triple[1]])
                data.append(triple[2])
            coboundary = coo_matrix((data, (row,col)), shape=(len(indices_of_actual_simplices[dim+1]),
                                                              len(indices_of_actual_simplices[dim])
                                                                )).toarray()
            coboundaries.append(coboundary)

        return coboundaries
    
class PersistentSheafLaplacian(Complex):
    """ Persistent Laplacian made with a cellular sheaf. The information is encoded in the sheaf_simplex_tree argument.
    
    """
    def __init__(self, sst: sheaf_simplex_tree, eigs_Algorithm = eigs_Algorithms.selfadjoint, up_Algorithm = up_Algorithms.schur):
        coboundaries, filtrations = sst.apply_restriction_function()
        boundaries = [x.T for x in coboundaries]    
        super().__init__(boundaries, filtrations, eigs_Algorithm, up_Algorithm,  storage = storage.float)


class NewPersistentSheafLaplacian(object):
    """ Persistent Laplacian made with a cellular sheaf. The information is encoded in the sheaf_simplex_tree argument.
    
    
    Attributes
    ----------
    verbose : boolean
        Print progress if spectra() is called
    flipped : boolean
        Compute the top-dimensional Laplacian's eigenvalues via the eigenvalues of the smaller of B_N B_N^T or B_N^T B_N and possible zero-padding
        
    Methods
    -------
    set_boundaries_filtrations(boundaries, filtrations)
        If the boundaries and filtrations were not set in the constructor, set them here.
    get_L(dim, a, b)
        Get the persistent Laplacian matrix.
    get_up(dim, a, b)
        Get the persistent up-Laplacian matrix.
    get_down(dim, a)
        Get the persistent down-Laplacian matrix.
    nonzero_spectra(dim, a, b, PH_basis=None, use_dummy_harmonic_basis=True)
        Compute the nonzero eigenvalues of the PL using reduction by persistent homology or standard basis of the null space.
    spectra(dim = None, a = None, b = None, request_list = None)
        Compute the eigenvalues of L_{dim}^{a,b} or for every tuple (dim, a, b) in request_list.
    eigenpairs(dim = None, a = None, b = None, request_list = None)
        Compute the eigenvalues and eigenvectors of L_{dim}^{a,b} or for every tuple (dim, a, b) in request_list.
    eigenvalues_summarize(eigenvalues)
        Compute the betti number and least nonzero eigenvalue of a list of eigenvalues.
    print_boundaries(self):
        Print all boundaries and corresponding filtrations in the complex.
    store_spectra(spectra_list, file_prefix):
        Store all of the eigenvalues in files f"{file_prefix}_spectra_{dim}.txt" for each dimension in the complex.    
    store_spectra_summary(spectra_list, file_prefix):
        Store eigenvalue summaries in file f"{file_prefix}_spectra_summary.txt".
    filtration_list_to_spectra_request(filtrations, dims):
        Get a list of tuples (dim, a, b) for all combinations of dimension and successive filtration values (a=filtrations[i], b=filtrations[i+1]).
    get_all_filtrations(self):
        Get a sorted list of all filtration values that occur in the complex.
    """

    # def __init__(self, sst: sheaf_simplex_tree, eigs_Algorithm = eigs_Algorithms.selfadjoint, up_Algorithm = up_Algorithms.schur):
        # coboundaries, filtrations = sst.apply_restriction_function()
        # boundaries = [x.T for x in coboundaries]    
        # super().__init__(boundaries, filtrations, eigs_Algorithm, up_Algorithm,  storage = storage.float)

    def __init__(self,sst: sheaf_simplex_tree = None, boundaries = None, filtrations = None, eigs_Algorithm = eigs_Algorithms.selfadjoint, up_Algorithm = up_Algorithms.schur, simplex_tree = None, storage = None):
        pl_class = getattr(_petls, "NewPersistentSheafLaplacian")

        if sst is None and boundaries is None and filtrations is None:
            raise TypeError("NewPersistentSheafLaplacian requires either (a) a sheaf_simplex_tree or (b) boundaries and filtrations. All were None.")

        if sst is not None:
            coboundaries, filtrations = sst.apply_restriction_function()
            boundaries = [x.T for x in coboundaries]

        self.pl = pl_class(boundaries,filtrations)
        self.verbose = False
        self.flipped = False


    def set_boundaries_filtrations(self, boundaries, filtrations):
        """ If the boundaries and filtrations were not set in the constructor, set them here.

        Parameters
        ----------
        boundaries : List[np.array]
            List of boundary matrices
        filtrations : List[List[float]]
            For each dimension in the complex, a list of simplex filtration values

        Returns
        -------
        None

        Examples
        --------
        >>> d1 = np.array([[-1,0,-1],
                [1,-1,0],
                [0,1,1]])
        >>> d2 = np.array([[1],[1],[-1]]) 
        >>> boundaries = [d1,d2]
        >>> filtrations = [[0,1,2],[3,4,5],[5]]
        >>> pl = petls.Complex()
        >>> pl.set_boundaries_filtrations(boundaries, filtrations)

        """
        self.pl.set_boundaries_filtrations(boundaries, filtrations)
    
    @property
    def verbose(self):
        return self._verbose

    @verbose.setter
    def verbose(self, verbose):
        self._verbose = verbose
        self.pl.set_verbose(verbose)

    @property
    def flipped(self):
        return self._flipped
    
    @flipped.setter
    def flipped(self, flipped):
        self._flipped = flipped
        self.pl.set_flipped(flipped)

    def set_flipped(self, flipped=True):
        self.pl.set_flipped(flipped)

    def get_L(self, dim, a, b):
        return self.pl.get_L(dim, a, b)
    
    def get_up(self, dim, a, b):
        return self.pl.get_up(dim, a, b)

    def get_down(self, dim, a):
        return self.pl.get_down(dim, a)

    def nonzero_spectra(self, dim, a, b, PH_basis=None, use_dummy_harmonic_basis=True):
        self.pl.nonzero_spectra(dim, a, b, PH_basis, use_dummy_harmonic_basis)


    def spectra(self,dim = None, a = None, b = None, request_list = None, allpairs = False):
        """ Compute the eigenvalues of L_{dim}^{a,b} or for every tuple (dim, a, b) in request_list.

        Parameters
        ----------
        dim : int, optional
            Dimension.
        a : float, optional
            Start filtration value.
        b : float, optional
            End filtration value
        request_list : List[List[int, float, float]], optional
            List of (dim, a, b) to compute the spectra of

        Returns
        -------
        List[float]
            If passed dim, a, and b, returns eigenvalues of L_{dim}^{a,b}
        List[Tuple[int, float, float, List[float]]]
            If passed request_list, returns a tuple (dim, a, b, eigenvalues) for each request in request_list    
            If passed no arguments, returns as if request_list is of all combinations (dim, a, b) where b is the next filtration after a, and (dim, a, a) when a is the largest filtration value  
        
        Examples
        --------
        >>> d1 = np.array([[-1,0,-1],
                [1,-1,0],
                [0,1,1]])
        >>> d2 = np.array([[1],[1],[-1]]) 
        >>> boundaries = [d1,d2]
        >>> filtrations = [[0,1,2],[3,4,5],[5]]
        >>> pl = Complex(boundaries, filtrations)
        >>> pl.spectra(0, 1.2, 4.5)
        [0.0,  1.9999998807907104]
        >>> pl.spectra(request_list = [[0, 1.2, 4.5]])
        [(0, 1.2, 4.5, [0.0, 1.9999998807907104])]
        >>> pl.spectra()
        [(0, 0.0, 3.0, [0.0]), (1, 0.0, 3.0, []), (2, 0.0, 3.0, []), (0, 3.0, 4.0, [0.0, 0.9999998807907104, 3.0]), (1, 3.0, 4.0, [2.0]), (2, 3.0, 4.0, []), (0, 4.0, 5.0, [0.0, 2.999999761581421, 3.0]), (1, 4.0, 5.0, [0.9999999403953552, 2.999999761581421]), (2, 4.0, 5.0, []), (0, 5.0, 5.0, [0.0, 2.999999761581421, 3.0]), (1, 5.0, 5.0, [3.0, 3.0, 3.0]), (2, 5.0, 5.0, [3.0])]
        """
        
        
        if request_list is not None:
            return self.pl.spectra(request_list)
        elif dim is not None and a is not None and b is not None:
            return self.pl.spectra(dim, a, b)
        elif allpairs:
            return self.pl.spectra_allpairs()
        else:
            return self.pl.spectra()
        
    def eigenpairs(self,dim = None, a = None, b = None, request_list = None):
        """ Compute the eigenvalues and eigenvectors of L_{dim}^{a,b} or for every tuple (dim, a, b) in request_list.

        Parameters
        ----------
        dim : int, optional
            Dimension.
        a : float, optional
            Start filtration value.
        b : float, optional
            End filtration value
        request_list : List[List[int, float, float]], optional
            List of (dim, a, b) to compute the spectra of

        Returns
        -------
        Tuple[List[float], numpy.ndarray]
            If passed dim, a, and b, returns eigenvalues and eigenvectors of L_{dim}^{a,b}
        List[Tuple[int, float, float, List[float], numpy.ndarray]]
            If passed request_list, returns a tuple (dim, a, b, eigenvalues, eigenvectors) for each request in request_list.   
            If passed no arguments, returns as if request_list is of all combinations (dim, a, b) where b is the next filtration after a, and (dim, a, a) when a is the largest filtration value.  
            If passed only one argument, assume it is a request_list.
        
        Examples
        --------
        >>> d1 = np.array([[-1,0,-1],
                [1,-1,0],
                [0,1,1]])
        >>> d2 = np.array([[1],[1],[-1]]) 
        >>> boundaries = [d1,d2]
        >>> filtrations = [[0,1,2],[3,4,5],[5]]
        >>> pl = Complex(boundaries, filtrations)
        >>> pl.spectra(0, 1.2, 4.5)
        [0.0,  1.9999998807907104]
        >>> pl.spectra(request_list = [[0, 1.2, 4.5]])
        [(0, 1.2, 4.5, [0.0, 1.9999998807907104])]
        >>> pl.spectra()
        [(0, 0.0, 3.0, [0.0]), (1, 0.0, 3.0, []), (2, 0.0, 3.0, []), (0, 3.0, 4.0, [0.0, 0.9999998807907104, 3.0]), (1, 3.0, 4.0, [2.0]), (2, 3.0, 4.0, []), (0, 4.0, 5.0, [0.0, 2.999999761581421, 3.0]), (1, 4.0, 5.0, [0.9999999403953552, 2.999999761581421]), (2, 4.0, 5.0, []), (0, 5.0, 5.0, [0.0, 2.999999761581421, 3.0]), (1, 5.0, 5.0, [3.0, 3.0, 3.0]), (2, 5.0, 5.0, [3.0])]
        """
        
        if (dim is not None and a is None and b is None and request_list is None):
            request_list = dim
            dim = None
        
        if request_list is not None:
            return self.pl.eigenpairs(request_list)
        elif dim is not None and a is not None and b is not None:
            return self.pl.eigenpairs(dim, a, b)
        else:
            return self.pl.eigenpairs()
        
    def eigenvalues_summarize(self, eigenvalues):
        return self.pl.eigenvalues_summarize(eigenvalues)
    
    def store_L(self, dim, a, b, prefix):
        # caution: recomputes L
        self.pl.store_L(dim, a, b, prefix)
    
    def print_boundaries(self):
        self.pl.print_boundaries()
    
    def store_spectra(self, spectra_list, file_prefix):
        self.pl.store_spectra(spectra_list, file_prefix)
    
    def store_spectra_summary(self, spectra_list, file_prefix):
        self.pl.store_spectra_summary(spectra_list, file_prefix)

    def filtration_list_to_spectra_request(self, filtrations, dims):
        return self.pl.filtration_list_to_spectra_request(filtrations, dims)
    
    def get_all_filtrations(self):
        return self.pl.get_all_filtrations()
    