# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'petls'
copyright = '2025, Benjamin Jones'
author = 'Benjamin Jones'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.autodoc',
    'sphinx.ext.coverage',
    # 'sphinx.ext.githubpages', # TODO: enable this once repository is public
    'sphinx.ext.viewcode',
    'sphinx.ext.mathjax',
    'numpydoc',
    'sphinx_copybutton',
    'sphinx_design',
    'breathe',
]

import sys
sys.path.append("/home/jones657/petls/docs/breathe/")
breathe_projects = {"petls": "/home/jones657/petls/docs/xml/"}
breathe_default_project = "petls"

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'alabaster'
html_theme = 'pydata_sphinx_theme'
# html_static_path = ['_static']
