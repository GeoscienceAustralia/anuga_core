# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------


# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here.
#import os
#import sys
#sys.path.insert(0, pathlib.Path(__file__).parents[2].resolve().as_posix())
#sys.path.insert(0, os.path.abspath('../..'))

# -- Project information -----------------------------------------------------

project = 'ANUGA'
copyright = 'Commonwealth of Australia (Geoscience Australia) and the Australian National University 2004-Now'
author =  'Stephen Roberts, Ole Nielsen, Gareth Davies'

# The full version, including alpha/beta/rc tags
import anuga
release = anuga.__version__


import os
import sys
sys.path.insert(0, os.path.abspath("../../anuga"))
sys.path.insert(0, os.path.abspath("_ext"))

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

import sphinx_rtd_theme

extensions = [
    'sphinx_rtd_theme',
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    'sphinx.ext.mathjax',
    'sphinx.ext.coverage',
    'sphinx.ext.viewcode',
    'sphinx.ext.autosectionlabel',
    'sphinx_copybutton',
    'nbsphinx',
    'spec_label',
    'myst_parser',      # renders the Markdown physics specification
]

# MathJax 4 breaks inline maths at relations and operators, so a paragraph can
# end with "S_s" and start the next line with "= E_s - D_s". Keep each inline
# expression on one line; long displayed equations may still break.
mathjax4_config = {'output': {'linebreaks': {'inline': False}}}

def linkcode_resolve(domain, info):
    if domain != 'py':
        return None
    if not info['module']:
        return None
    filename = info['module'].replace('.', '/')
    return "https://somesite/sourcerepo/%s.py" % filename


#autodoc_mock_imports = ["anuga"]

# Document class members by default so the method/attribute summary tables on
# each class page link through to the individual method signatures + docstrings.
# (autodoc_default_flags is deprecated and ignored by modern Sphinx — use
# autodoc_default_options.)
autodoc_default_options = {
    'members': True,
    'show-inheritance': True,
}
autosummary_generate = True

# Render NumPy-style 'Attributes' sections as :ivar: fields rather than
# .. attribute:: directives, so an attribute that is also a property (e.g.
# Geo_reference.epsg) is not documented twice (avoids a duplicate-object warning).
napoleon_use_ivar = True

autosectionlabel_prefix_document = True

# Suppress duplicate-label warnings from autosectionlabel on autodoc-generated
# docstring sections (e.g. "Parameters", "Returns" in multiple functions on
# the same page).  Explicit .. _label: targets are used for all real
# cross-references so these auto-labels are not needed.
suppress_warnings = ['autosectionlabel']

#extensions.append('sphinxcontrib.bibtex')
#bibtex_bibfiles = ['refs.bib']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Keeps the bracketed physics labels from breaking mid-token; see the file.
html_css_files = ['anuga-spec-labels.css']