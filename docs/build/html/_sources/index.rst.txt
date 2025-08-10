.. PersistentLaplacians documentation master file, created by
   sphinx-quickstart on Wed Jun 12 17:44:33 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PersistentLaplacians's documentation!
================================================


.. note::
   This project is under active development as of January 2025. For developing it, I have restricted some of the configurations for usage and installation to my environment. If you want a minor change to setup (e.g. build wheels for Python 3.11), just ask!

Overview
********

Persistent Laplacians are a persistent version of combinatorial Laplacians, which generalize the well-studied graph Laplacian to simplicial and other complexes, such as path complexes. As with graph Laplacians, we study persistent Laplacians largely through their eigenvalues.

This is a **C++ library** with **Python bindings** that computes Persistent Laplacians extremely quickly (orders of magnitude faster than other tools) for a variety of simplicial and non-simplicial complexes. Most importantly, it will compute the persistent Laplacian given any collection of (filtered) boundary matrices. 

This project is intended to be user-friendly for two audiences:

1. Researchers in Topological Data Analysis (TDA), who develop and analyze both theoretical and computational tools.
2. Researchers with data that could be analyzed with topology. If you use Persistent Homology, it is likely that Persistent Laplacians will be of interest to you!

Source code is available on `GitHub <https://github.com/bdjones13/PersistentLaplacians/>`_.

This project is written and maintained by `Ben Jones <https://www.benjones-math.com/>`_ at Michigan State University. 




.. grid:: 2

    .. grid-item-card::

        Python Documentation
        ^^^^^^^^^^^^^^^^^^^^

        Python 

        +++

        .. button-ref:: python
            :expand:
            :color: secondary
            :click-parent:

            To the python docs

    .. grid-item-card::

        C++ Documentation
        ^^^^^^^^^^^^^^^^^

        Documentation for C++

        +++

        .. button-ref:: cpp
            :expand:
            :color: secondary
            :click-parent:

            To the C++ docs

Acknowledgements:
*****************

This project is partially supported by grants xxx.

.. toctree::
   :maxdepth: 1
   :hidden:
   :titlesonly:

   python/index
   cpp/index
   Getting_Started
