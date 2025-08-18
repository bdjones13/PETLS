See [the documentation](https://benjones-math.com/software/PETLS/) for details on requirements, installation, and API usage.

# Python Dependencies

* Eigen version 3.4 or higher
* pybind11
* scipy
* numpy
* Boost

# Python installation and usage
The simplest way to install is:

``` pip install petls```

The next simplest is to clone the repository and run

``` pip install .```


After that, you're all set:
```
   $ python
   >>> from petls import Rips
   >>> import numpy as np
   >>> points = np.array([[0, 0], [0, 3], [4, 0], [4, 3]])
   >>> complex = Rips(points, max_dim=3)
   >>> complex.spectra(dim=1, a=3, b=4)
   [2.0, 2.0]

```

# C++ installation
First, make and enter the ```cpp/build``` directory.

```
   cd cpp
   mkdir build
   cd build
```

Next, we make and install the project. If you have root access:
```
   cmake ..
   sudo make all install
```

If you do not have root access (e.g. on an HPC):
```
   cmake -DCMAKE_INSTALL_PREFIX=<location> ..
   make all install
```

Usually ```CMAKE_INSTALL_PREFIX``` is something like ```~/custom_libraries/```. 

For examples of how to include the library, see the ```examples``` directory, e.g. ```examples/basic_cpp/```.

# Acknowledgements
This work was supported in part by NIH grant R35GM148196, NSF grant DMS-2052983, and MSU Research Foundation.

# Citation
If you wish to cite this work, please use the following citation:
```
@misc{jones2025petlspersistenttopologicallaplacian,
      title={PETLS: PErsistent Topological Laplacian Software}, 
      author={Benjamin Jones and Guo-Wei Wei},
      year={2025},
      eprint={2508.11560},
      archivePrefix={arXiv},
      primaryClass={math.AT},
      url={https://arxiv.org/abs/2508.11560}, 
}
```