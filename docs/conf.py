# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import html
import os
import re
import sys

# Add the project root to the path for autodoc
sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------

project = 'mlmm'
copyright = '2025, Takuto Ohmura'
author = 'Takuto Ohmura'

# Hardcoded release version (setuptools_scm is unavailable in the docs build env)
release = '0.4.0'

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

# Keep "--flag" in link text as typed (no en/em-dash conversion).
smartquotes_action = 'qe'

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

# Same palette and page design as the pdb2reaction docs (the p2r-* names of the CSS
# classes and of these settings are shared with them so the two stay in step): navy
# #2E5B96, light navy #7BA3D7, stage colours (MEP #70AD47, TS #4747C1, IRC #2E9691,
# thermo #EC7C30), and a rose for DFT.
# In the light theme green/teal/orange are darkened so text in them reaches 4.5:1.
_FONTS = {
    'font-stack': '"Inter", "Noto Sans JP", -apple-system, BlinkMacSystemFont, "Segoe UI", "Hiragino Sans", "Yu Gothic UI", Meiryo, sans-serif',
    'font-stack--headings': '"Inter", "Noto Sans JP", -apple-system, BlinkMacSystemFont, "Segoe UI", "Hiragino Sans", "Yu Gothic UI", Meiryo, sans-serif',
    'font-stack--monospace': '"JetBrains Mono", ui-monospace, SFMono-Regular, Menlo, Consolas, monospace',
}

# Reaction-workflow stages shown as the stage strip of the top page and the stage
# badges of the command pages. Labels follow the steps of "How it works" in
# all.md / ja/all.md: step 3 is split into scan and path, step 4 into TS and IRC, and
# step 5 into thermochemistry and DFT (one icon each, _static/icons/<p2r_icon_set>/);
# each stage takes its colour from the palette above.
# pages[0] is the page a stage links to; commands not listed here stay neutral.
p2r_pipeline = [
    # Steps 1 and 2 of "How it works" (ML region, then topology and layers) form one stage.
    {'id': 'prep', 'en': 'ML/MM model setup', 'ja': 'ML/MM モデルのセットアップ',
     'pages': ['extract', 'mm-parm', 'define-layer', 'fix-altloc', 'add-elem-info'], 'color': 'navy', 'icon': 'extract'},
    # Scan comes right before the MEP search; only the Scan-list mode goes through it.
    {'id': 'scan', 'en': 'Scan', 'ja': 'スキャン',
     'pages': ['scan', 'scan2d', 'scan3d'], 'color': 'scan', 'icon': 'scan'},
    {'id': 'path', 'en': 'Building the path', 'ja': '経路の作成',
     'pages': ['path-opt', 'path-search'], 'color': 'green', 'icon': 'mep'},
    {'id': 'ts', 'en': 'Optimizing the TS', 'ja': 'TS の最適化',
     'pages': ['tsopt'], 'color': 'violet', 'icon': 'tsopt'},
    {'id': 'irc', 'en': 'Following the IRC', 'ja': 'IRC の追跡',
     'pages': ['irc'], 'color': 'teal', 'icon': 'irc'},
    # \u00ad (soft hyphen) lets the narrow strip column break it as "Thermo-chemistry".
    {'id': 'thermo', 'en': 'Thermo\u00adchemistry', 'ja': '熱化学',
     'pages': ['freq'], 'color': 'orange', 'icon': 'thermo'},
    {'id': 'dft', 'en': 'DFT', 'ja': 'DFT',
     'pages': ['dft'], 'color': 'rose', 'icon': 'dft'},
]
# Stage icon set: the folder under _static/icons/ (u50 = every line drawn 5.0 units
# wide on the 96-unit canvas).
p2r_icon_set = 'u50'
# Stages of each numbered step of "How it works" in all.md (steps 1 and 2 share one stage;
# steps 3-5 cover two each).
p2r_all_steps = [['prep'], ['prep'], ['scan', 'path'], ['ts', 'irc'], ['thermo', 'dft']]
p2r_pipeline_all = 'all'      # the page that runs every stage (label: its title)
p2r_pipeline_all_label = {'en': 'End-to-end workflow', 'ja': '一気通貫ワークフロー'}

# Stage decorations: stage dots in the sidebar, tables and headings, and the hairline in the
# stage colours under the top bar (pipeline.css, under the html class p2r-deco).
p2r_deco = True

# Stages each input mode of `all` goes through, shown in the entry cards of the top
# page (all.md "How it works"; quickstart-*.md; TS-only mode skips the MEP search).
p2r_mode_stages = {
    'endpoint': ['prep', 'path', 'ts', 'irc', 'thermo', 'dft'],
    'scan': ['prep', 'scan', 'path', 'ts', 'irc', 'thermo', 'dft'],
    'tsonly': ['prep', 'ts', 'irc', 'thermo', 'dft'],
}

html_context = {
    'p2r_repo_url': 'https://github.com/t-0hmura/mlmm_toolkit',   # GitHub icon of the top bar
    'p2r_deco': p2r_deco,
    'p2r_mode_stages': p2r_mode_stages,
    'p2r_pipeline': p2r_pipeline,
    'p2r_pipeline_all': p2r_pipeline_all,
    'p2r_pipeline_all_label': p2r_pipeline_all_label,
    'p2r_all_steps': p2r_all_steps,
    'p2r_stage_of': {**{p: s['id'] for s in p2r_pipeline for p in s['pages']},
                     p2r_pipeline_all: 'all'},
}

_STAGE_VARS = {
    **{f"p2r-{s['id']}": f"var(--p2r-{s['color']})" for s in p2r_pipeline},
    'p2r-gradient': 'linear-gradient(90deg, %s)' % ', '.join(
        f"var(--p2r-{s['id']})" for s in p2r_pipeline),
}

html_theme_options = {
    'light_css_variables': {
        **_FONTS,
        **_STAGE_VARS,
        'p2r-navy': '#2E5B96',
        'p2r-navy-deep': '#1B365D',
        'p2r-sky': '#7BA3D7',
        'p2r-violet': '#4444C0',
        'p2r-teal': '#217A76',
        'p2r-green': '#4A7F2A',
        'p2r-scan': '#566A86',   # slate (scan colour S5)
        'p2r-orange': '#B4520F',
        'p2r-rose': '#B0396E',
        'p2r-on-stage': '#FFFFFF',
        'p2r-tint': '#EEF3FA',
        'p2r-blue': '#1F6FEB',   # the dot of `all` in the sidebar and tables
        'p2r-tint-strong': '#DCE7F5',
        'p2r-hero-from': '#F3F7FD',
        'p2r-hero-to': '#E3ECF8',
        'p2r-code-border': '#0F1D33',
        'p2r-btn-from': '#3E6CA8',
        'p2r-btn-to': '#1B365D',
        'p2r-figure-filter': 'none',
        'p2r-shadow': '0 1px 2px rgba(27, 54, 93, 0.06), 0 4px 14px rgba(27, 54, 93, 0.07)',
        'p2r-topbar-bg': 'rgba(255, 255, 255, 0.82)',
        'p2r-surface': 'rgba(255, 255, 255, 0.88)',
        'p2r-border-strong': '#C9D5E6',
        'color-brand-primary': '#2E5B96',
        'color-brand-content': '#2E5B96',
        'color-brand-visited': '#3D4FA8',
        'color-foreground-primary': '#18222F',
        'color-foreground-secondary': '#4A5568',
        'color-foreground-muted': '#68748A',
        'color-foreground-border': '#8A96AA',
        'color-background-primary': '#FFFFFF',
        'color-background-secondary': '#F5F8FC',
        'color-background-hover': '#E8EFF8',
        'color-background-hover--transparent': '#E8EFF800',
        'color-background-border': '#DDE5F0',
        'color-background-item': '#C7D3E3',
        'color-link': '#2E5B96',
        'color-link--hover': '#4444C0',
        'color-link-underline': 'rgba(46, 91, 150, 0.28)',
        'color-link-underline--hover': '#4444C0',
        'color-sidebar-background': '#F5F8FC',
        'color-sidebar-background-border': '#DDE5F0',
        'color-sidebar-brand-text': '#1B365D',
        'color-sidebar-caption-text': '#2E5B96',
        'color-sidebar-link-text': '#334155',
        'color-sidebar-link-text--top-level': '#1F2D40',
        'color-sidebar-item-background--current': '#E3ECF8',
        'color-sidebar-item-background--hover': '#E8EFF8',
        'color-sidebar-search-background': '#FFFFFF',
        'color-sidebar-search-border': '#DDE5F0',
        'color-toc-item-text--active': '#2E5B96',
        'color-table-header-background': '#EEF3FA',
        'color-table-border': '#DDE5F0',
        'color-inline-code-background': '#EEF3FA',
        'color-code-background': '#0F1D33',
        'color-code-foreground': '#E4EBF5',
        'color-highlight-on-target': '#FFF6D6',
        'color-admonition-background': '#FFFFFF',
        'color-admonition-title--note': '#2E5B96',
        'color-admonition-title-background--note': '#EEF3FA',
        'color-admonition-title--tip': '#217A76',
        'color-admonition-title-background--tip': '#E8F5F4',
        'color-admonition-title--important': '#4444C0',
        'color-admonition-title-background--important': '#EEEEFA',
        'color-admonition-title--warning': '#B4520F',
        'color-admonition-title-background--warning': '#FDF0E6',
    },
    'dark_css_variables': {
        **_FONTS,
        **_STAGE_VARS,
        'p2r-navy': '#8DB3E6',
        'p2r-navy-deep': '#DCE7F7',
        'p2r-sky': '#7BA3D7',
        'p2r-violet': '#A3A3F0',
        'p2r-teal': '#5CC4BE',
        'p2r-green': '#98CB6F',
        'p2r-scan': '#A9B8CE',   # slate (scan colour S5)
        'p2r-orange': '#F2A05E',
        'p2r-rose': '#F29BC2',
        'p2r-on-stage': '#0D131F',
        'p2r-tint': '#16233A',
        'p2r-blue': '#58A6FF',
        'p2r-tint-strong': '#1D2D49',
        'p2r-hero-from': '#121D31',
        'p2r-hero-to': '#0F182A',
        'p2r-code-border': '#25344F',
        'p2r-btn-from': '#5585C8',
        'p2r-btn-to': '#2E5B96',
        'p2r-figure-filter': 'brightness(0.93)',
        'p2r-shadow': '0 1px 2px rgba(0, 0, 0, 0.3), 0 6px 18px rgba(0, 0, 0, 0.25)',
        'p2r-topbar-bg': 'rgba(13, 19, 31, 0.82)',
        'p2r-surface': 'rgba(22, 35, 58, 0.72)',
        'p2r-border-strong': '#2C3D5C',
        'color-brand-primary': '#8DB3E6',
        'color-brand-content': '#8DB3E6',
        'color-brand-visited': '#B0B0F2',
        'color-foreground-primary': '#E3E9F2',
        'color-foreground-secondary': '#AAB6C9',
        'color-foreground-muted': '#8592A8',
        'color-foreground-border': '#5D6A80',
        'color-background-primary': '#0D131F',
        'color-background-secondary': '#111A29',
        'color-background-hover': '#18243A',
        'color-background-hover--transparent': '#18243A00',
        'color-background-border': '#223049',
        'color-background-item': '#3A4A66',
        'color-link': '#8DB3E6',
        'color-link--hover': '#B0B0F2',
        'color-link-underline': 'rgba(141, 179, 230, 0.32)',
        'color-link-underline--hover': '#B0B0F2',
        'color-sidebar-background': '#101827',
        'color-sidebar-background-border': '#1E2B42',
        'color-sidebar-brand-text': '#E3E9F2',
        'color-sidebar-caption-text': '#8DB3E6',
        'color-sidebar-link-text': '#B9C4D6',
        'color-sidebar-link-text--top-level': '#D5DDEA',
        'color-sidebar-item-background--current': '#18263F',
        'color-sidebar-item-background--hover': '#18243A',
        'color-sidebar-search-background': '#0D131F',
        'color-sidebar-search-border': '#223049',
        'color-toc-item-text--active': '#8DB3E6',
        'color-table-header-background': '#16233A',
        'color-table-border': '#223049',
        'color-inline-code-background': '#18263F',
        'color-code-background': '#0A1222',
        'color-code-foreground': '#E4EBF5',
        'color-highlight-on-target': '#2E2A14',
        'color-admonition-background': '#111A29',
        'color-admonition-title--note': '#8DB3E6',
        'color-admonition-title-background--note': '#16233A',
        'color-admonition-title--tip': '#5CC4BE',
        'color-admonition-title-background--tip': '#122B2E',
        'color-admonition-title--important': '#A3A3F0',
        'color-admonition-title-background--important': '#1C1E3D',
        'color-admonition-title--warning': '#F2A05E',
        'color-admonition-title-background--warning': '#2E2117',
    },
    'sidebar_hide_name': False,
    'navigation_with_keys': True,
}

# Code blocks are dark navy in both themes.
pygments_style = 'github-dark'
pygments_dark_style = 'github-dark'

html_title = 'mlmm-toolkit'
html_short_title = 'mlmm-toolkit'

# Static files path
html_static_path = ['_static']

# Custom JS — override Furo 3-state toggle to 2-state (light ↔ dark)
html_js_files = ['theme_toggle.js', 'pipeline.js']

# Custom CSS (Inter / Noto Sans JP / JetBrains Mono from Google Fonts);
# p2r-icons.css is written at the end of the build by _p2r_icon_css below.
html_css_files = [
    'https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800'
    '&family=JetBrains+Mono:wght@400;500;600&family=Noto+Sans+JP:wght@400;500;700'
    '&display=swap',
    'custom.css',
    'pipeline.css',
    'p2r-icons.css',
]

# Favicon: the mlmm-toolkit icon (a placeholder until the logo is designed)
html_favicon = '_static/mlmm-toolkit-icon.svg'

# Logo (optional)
# html_logo = '_static/logo.png'

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


# -- Stage icons ---------------------------------------------------------------

def _p2r_icon_css(app, exception):
    """Write _static/p2r-icons.css: each stage icon of _static/icons/<p2r_icon_set>/ as a CSS mask.

    The SVGs are embedded as data URIs because browsers fetch mask images in CORS
    mode, so a url() to a file fails when the pages are opened from file://.
    """
    if exception is not None or app.builder.format != 'html':
        return
    from pathlib import Path
    from urllib.parse import quote
    rules = [f'/* Generated by conf.py from _static/icons/{p2r_icon_set}/*.svg (stage icons as masks). */']
    for s in p2r_pipeline:
        svg = (Path(app.srcdir) / '_static' / 'icons' / p2r_icon_set / f"{s['icon']}.svg").read_text(encoding='utf-8')
        uri = 'data:image/svg+xml,' + quote(' '.join(svg.split()).replace('"', "'"), safe=" '=:/.,-;()")
        rules.append(f'.p2r-ico--{s["id"]} {{ --p2r-ico: url("{uri}"); }}')
    (Path(app.outdir) / '_static' / 'p2r-icons.css').write_text('\n'.join(rules) + '\n', encoding='utf-8')


# The <code> part stays inside one link, so a match never runs on into the next entries.
_NAV_CMD = re.compile(r'(<a class="[^"]*reference internal[^"]*" href="[^"]*")(>)(<code[^>]*>(?:(?!</a>).)*?</code>)\s*([（(][^<]*[）)])(</a>)', re.S)

# Hover text of the quick-start entries, whose mode names alone do not say what they are for
# (the purposes on the entry cards of the top page).
p2r_nav_purpose = {
    'quickstart-all': 'Analyze the mechanism end to end from the structures before and after the reaction',
    'quickstart-scan': 'Analyze the mechanism end to end from one structure',
    'quickstart-tsopt': 'Analyze the mechanism end to end from a TS structure',
    'ja/quickstart-all': '反応の前後の構造から反応機構解析を一気通貫で行う',
    'ja/quickstart-scan': '1 つの構造から一気通貫で反応機構解析を行う',
    'ja/quickstart-tsopt': 'TS 構造から一気通貫で反応機構解析を行う',
}


def _p2r_sidebar(app, pagename, templatename, context, doctree):
    """Sidebar: a Home item on top, command entries shown as the command name only, and the
    purpose of each quick-start entry as its hover text (p2r_nav_purpose).

    A command entry "`extract`（...）" keeps its full title as the hover text (title attribute);
    the description stays in the markup but is hidden by CSS (.p2r-nav-desc).
    """
    tree = context.get('furo_navigation_tree')
    if not tree:
        return

    def short(m):
        full = html.unescape(re.sub(r'<[^>]+>', '', m.group(3)) + m.group(4))
        return (f'{m.group(1)} title="{html.escape(full, quote=True)}"{m.group(2)}{m.group(3)}'
                f'<span class="p2r-nav-desc">{m.group(4)}</span>{m.group(5)}')

    tree = _NAV_CMD.sub(short, tree)
    for doc, text in p2r_nav_purpose.items():
        href = '#' if doc == pagename else context['pathto'](doc)
        tree = tree.replace(f'href="{href}">', f'href="{href}" title="{html.escape(text, quote=True)}">', 1)
    is_ja = pagename.startswith('ja/')
    home = 'ja/index' if is_ja else master_doc
    cur = ' current current-page' if pagename == home else ''
    label = 'ホーム' if is_ja else 'Home'
    href = context['pathto'](home)
    context['furo_navigation_tree'] = (
        f'<ul class="p2r-home-list"><li class="toctree-l1 p2r-home{cur}">'
        f'<a class="reference internal{" current" if cur else ""}" href="{href}">{label}</a></li></ul>' + tree)


def setup(app):
    app.connect('build-finished', _p2r_icon_css)
    # After Furo's own html-page-context handler, which builds furo_navigation_tree.
    app.connect('html-page-context', _p2r_sidebar, priority=900)
