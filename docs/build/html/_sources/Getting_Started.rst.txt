
Getting Started
=======================================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. Note::
   Using Intel MKL or BLAS offers a speedup by a factor of about 6. This documentation does not currently reflect how to do that.

Python Installation
*******************
If you have a system compatible with any of the binary wheels listed `here <https://test.pypi.org/project/PersistentLaplacians/>`_, you can just install via ``pip install PersistentLaplacians``.

Otherwise, to install from source, there are the following dependencies:

- CMake >= 3.16.3 
- Python >= 3.10
- pytest

If you intend to use the Alpha complex from Gudhi, you will also need the following dependencies at the time of installation:

- Boost >= 1.78.0
- CGAL >= 4.11

There are three ways to install from source:

1. If you do not have a system compatible with any of the binaries on PyPI, then ``pip install PersistentLaplacians`` should do the job.
2. Clone the `GitHub repository <https://github.com/bdjones13/PersistentLaplacians/>`_ and from the project root run ``pip install .``
3. Clone the `GitHub repository <https://github.com/bdjones13/PersistentLaplacians/>`_ and from the project root run:: 

      mkdir build
      cd build
      cmake ..
      make
      sudo make install


C++ Installation
****************
Dependencies:

- CMake >= 3.16.3 
- Eigen 3 

.. note::
   
   This project downloads and compiles Eigen 3.4 for internal usage (some features new to 3.4 are used), but to call PersistentLaplacian functions you must pass in Eigen matrices, so you must have access to your own version of Eigen 3. You likely do not need Eigen 3.4 to use this library.

If you intend to use the Alpha complex from Gudhi, you will also need the following dependencies at the time of installation:

- Boost >= 1.78.0
- CGAL >= 4.11

There is one way to install from source. From the project root::

   cd cpp
   mkdir build
   cd build
   cmake ..
   make
   sudo make install


Usage
*****

Suppose you have the following filtered simplicial complex:

Dimension 0:

- point a, added at filtration = 0
- point b, added at filtration = 1
- point c, added at filtration = 2

Dimension 1:

- line (a,b), added at filtration = 3
- line (b,c), added at filtration = 4
- line (a,c), added at filtration = 5

Dimension 2:

- triangle (a,b,c), added at filtration = 5

You can create the persistent Laplacians and compute the spectra:

**Python**

.. code-block:: python

   import numpy as np
   import PersistentLaplacians


   # boundary matrices
   d1 = np.array([[-1,0,-1],
                [1,-1,0],
                [0,1,1]])
   d2 = np.array([[1],[1],[-1]]) 
   boundaries = [d1,d2]

   filtrations = [[0,1,2],   # dim 0 filtrations
                  [3,4,5],    # dim 1 filtrations
                  [5]]        # dim 2 filtrations

   pl = PersistentLaplacians.PersistentLaplacian(boundaries, filtrations)
   print(pl.spectra())


**C++**

.. code-block:: c

   #include "PersistentLaplacian.hpp"
   #include <Eigen/Dense>
   #include <vector>
   #include <iostream>

   SparseMatrixInt d1(3,3);
   SparseMatrixInt d2(3,1);
   d1.coeffRef(0,0) = -1;
   d1.coeffRef(0,1) = 0;
   d1.coeffRef(0,2) = -1;
   d1.coeffRef(1,0) = 1;
   d1.coeffRef(1,1) = -1;
   d1.coeffRef(1,2) = 0;
   d1.coeffRef(2,0) = 0;
   d1.coeffRef(2,1) = 1;
   d1.coeffRef(2,2) = 1;

   d2.coeffRef(0,0) = 1;
   d2.coeffRef(1,0) = 1;
   d2.coeffRef(2,0) = -1;

   std::vector<filtration_type> c0_filtrations = {0.0, 1.0, 2.0};
   std::vector<filtration_type> c1_filtrations = {3.0, 4.0, 5.0};
   std::vector<filtration_type> c2_filtrations = {5.0};

   std::vector<SparseMatrixInt> boundaries;
   boundaries.push_back(d1);
   boundaries.push_back(d2);

   std::vector<std::vector<filtration_type>> filtrations;
   filtrations.push_back(c0_filtrations);
   filtrations.push_back(c1_filtrations);
   filtrations.push_back(c2_filtrations);
   
   PersistentLaplacians::PersistentLaplacian pl(boundaries,filtrations);
   std::cout << pl.spectra() << std::endl;




