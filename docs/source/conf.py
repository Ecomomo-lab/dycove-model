# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'DYCOVE'
copyright = '2025, N. Tull, M. Brückner'
author = 'N. Tull, M. Brückner'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx.ext.graphviz',
    'sphinx.ext.mathjax',
    'sphinx.ext.githubpages',
    "matplotlib.sphinxext.plot_directive",
    'sphinx.ext.viewcode',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# autodoc_mock_imports = [
#     'anuga',
#     'bmi',
#     'bmi_python',
#     'netCDF4',
#     'imageio',
#     'tqdm',
#     'scipy',
#     'numpy',
#     'matplotlib',
# ]

# suppress_warnings = [
#     'autodoc.duplicate_object',
# ]

autosummary_generate = True
automodapi_inheritance_diagram = False
autodoc_default_options = {
    'members': True,
    'inherited-members': True,
    'private-members': False,
    #'undoc-members': False,
    #'show-inheritance': True,
    #'imported-members': False,
}

napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
