.. PETLS documentation master file, created by
   sphinx-quickstart on Wed Jun 12 17:44:33 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PETLS: PErsistent Topological Laplacian Software
================================================

Overview
********

Persistent Laplacians are a persistent version of combinatorial Laplacians, which generalize the well-studied graph Laplacian to simplicial and other complexes, such as path complexes. As with graph Laplacians, we study persistent Laplacians largely through their eigenvalues.

This is a **C++ library** with **Python bindings** that computes Persistent Laplacians extremely quickly (orders of magnitude faster than other tools) for a variety of simplicial and non-simplicial complexes. Most importantly, it will compute the persistent Laplacian given any collection of (filtered) boundary matrices. 

This project is intended for two audiences:

1. Researchers in Topological Data Analysis (TDA), who develop and analyze both theoretical and computational tools.
2. Researchers with data that could be analyzed with topology. If you use Persistent Homology, it is likely that Persistent Laplacians will be of interest to you!

Source code is available on `GitHub <https://github.com/bdjones13/PETLS/>`_.

This project is written and maintained by `Ben Jones <https://www.benjones-math.com/>`_ at Michigan State University. 




.. grid:: 2

    .. grid-item-card::

        .. button-ref:: python
            :expand:
            :color: secondary
            :click-parent:


            Click here for Python documentation

    .. grid-item-card::

        .. button-ref:: cpp
            :expand:
            :color: secondary
            :click-parent:


            Click here for C++ documentation     

Acknowledgements:
*****************

This work was supported in part by NIH grant R35GM148196, NSF grant DMS-2052983, and MSU Research Foundation.

.. toctree::
   :maxdepth: 1
   :hidden:
   :titlesonly:

   python/index
   cpp/index
   Getting_Started

Citation
********
Please use the following bibtex entry to cite this work:

.. code-block::

    @misc{jones2025petlspersistenttopologicallaplacian,
        title={PETLS: PErsistent Topological Laplacian Software}, 
        author={Benjamin Jones and Guo-Wei Wei},
        year={2025},
        eprint={2508.11560},
        archivePrefix={arXiv},
        primaryClass={math.AT},
        url={https://arxiv.org/abs/2508.11560}, 
    }
