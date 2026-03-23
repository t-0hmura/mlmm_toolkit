# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# Add the project root to the path for autodoc
sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------

project = 'mlmm'
copyright = '2025, Takuto Ohmura'
author = 'Takuto Ohmura'

# Get version from package (setuptools_scm)
try:
    from mlmm import __version__
    release = __version__
except ImportError:
    # Fallback for docs build without package installed
    try:
        from mlmm._version import __version__
        release = __version__
    except ImportError:
        release = '0.0.0.dev0'

version = release

# Extract short version (e.g., "0.1.2" from "0.1.2.dev136+gb070dbf49")
short_version = release.split('.dev')[0] if '.dev' in release else release.split('+')[0]

# -- General configuration ---------------------------------------------------

extensions = [
    'myst_parser',                    # Markdown support
    'sphinx.ext.autodoc',             # API documentation from docstrings
    'sphinx.ext.autosummary',         # Generate autodoc summaries
    'sphinx.ext.napoleon',            # Google/NumPy style docstrings
    'sphinx.ext.viewcode',            # Add links to source code
    'sphinx.ext.intersphinx',         # Link to other projects' docs
    'sphinx_copybutton',              # Copy button for code blocks
]

# MyST Parser configuration
myst_enable_extensions = [
    'colon_fence',      # ::: directives
    'deflist',          # Definition lists
    'html_admonition',  # HTML-style admonitions
    'html_image',       # HTML image syntax
    'substitution',     # Substitution syntax
    'tasklist',         # Task lists
    'attrs_inline',     # Inline {#id} attributes
    'attrs_block',      # Block {#id} attributes for headings
]

myst_heading_anchors = 3

# MyST substitutions for version display in Markdown files
# Use {{ version }} or {{ release }} in .md files
myst_substitutions = {
    'version': short_version,
    'release': release,
    'project': project,
}

# Source file suffixes
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

# The master toctree document
master_doc = 'index'

# Patterns to exclude
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'README.md']

# Templates path
templates_path = ['_templates']

# -- Options for HTML output -------------------------------------------------

html_theme = 'furo'

html_theme_options = {
    'light_css_variables': {
        'color-brand-primary': '#2980B9',
        'color-brand-content': '#2980B9',
    },
    'dark_css_variables': {
        'color-brand-primary': '#56B4E9',
        'color-brand-content': '#56B4E9',
    },
    'sidebar_hide_name': False,
    'navigation_with_keys': True,
}

html_title = 'mlmm-toolkit'
html_short_title = 'mlmm-toolkit'

# Static files path
html_static_path = ['_static']

# Custom JS — override Furo 3-state toggle to 2-state (light ↔ dark)
html_js_files = ['theme_toggle.js']

# Custom CSS
html_css_files = ['custom.css']

# -- Options for autodoc -----------------------------------------------------

autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'special-members': '__init__',
    'undoc-members': True,
    'exclude-members': '__weakref__',
}

autodoc_typehints = 'description'
autodoc_typehints_format = 'short'

# -- Options for intersphinx -------------------------------------------------

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'torch': ('https://pytorch.org/docs/stable/', None),
}

# -- Options for copy button -------------------------------------------------

copybutton_prompt_text = r'>>> |\.\.\. |\$ |> '
copybutton_prompt_is_regexp = True

# -- Napoleon settings -------------------------------------------------------

napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = False
napoleon_type_aliases = None
napoleon_attr_annotations = True
