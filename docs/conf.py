"""Sphinx configuration for the shachen documentation.

Build with ``uv run --group docs sphinx-build -b html docs docs/_build/html``
(or ``make -C docs html``). The version is read from pyproject.toml so it
never drifts from the package.
"""

import os
import tomllib
from pathlib import Path

_PYPROJECT = tomllib.loads((Path(__file__).resolve().parents[1] / "pyproject.toml").read_text())

project = "shachen"
author = "Han Xiao"
copyright = "2026, Han Xiao — Apache-2.0"

# The number in pyproject.toml is the truth for a tagged build: that build
# really is that release. It is a half-truth for the site published from main,
# which is usually a few commits ahead of the last tag — a reader comparing the
# page against the version they installed from PyPI would be reading about code
# that has not shipped. So the docs workflow passes a PEP 440 local version,
# `0.2.1+g3fee821`, naming the commit the site was actually built from; nothing
# else (a plain `make -C docs html`, a PR build) sets it, and the bare version
# stands.
version = _PYPROJECT["project"]["version"]
release = os.environ.get("SHACHEN_DOCS_RELEASE") or version

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "myst_parser",
    "sphinx_copybutton",
]

exclude_patterns = ["_build", ".DS_Store"]

# -- autodoc -----------------------------------------------------------------

# Source order matches the order the equations run in, which is how every
# module is written; alphabetical would scramble the algorithm.
autodoc_member_order = "bysource"
autodoc_default_options = {
    "members": True,
    "show-inheritance": True,
    "member-order": "bysource",
}
# Type hints stay in the signature: the docstrings describe arguments in prose
# rather than in :param: fields, so a generated "Parameters" block would list
# types with no text next to them.
autodoc_typehints = "signature"
# `#:` comments on the dataclass fields in constants.py carry the paper's
# equation numbers; without this autodoc drops the undocumented siblings.
autodoc_class_signature = "separated"
# satpy is the only module-level import of an optional extra (io.satellite);
# cartopy, matplotlib, earthaccess and s3fs are imported inside functions and
# so need no mock.
autodoc_mock_imports = ["satpy"]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "scipy": ("https://docs.scipy.org/doc/scipy", None),
    "xarray": ("https://docs.xarray.dev/en/stable", None),
    "pyresample": ("https://pyresample.readthedocs.io/en/stable", None),
}

# -- MyST --------------------------------------------------------------------

myst_enable_extensions = ["colon_fence", "deflist", "dollarmath"]
myst_heading_anchors = 3

# -- HTML --------------------------------------------------------------------

html_theme = "furo"
html_title = f"shachen {release}"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_theme_options = {
    # Dust yellow, the colour the algorithm modulates into the imagery.
    "light_css_variables": {
        "color-brand-primary": "#8a6100",
        "color-brand-content": "#8a6100",
    },
    "dark_css_variables": {
        "color-brand-primary": "#e8b83c",
        "color-brand-content": "#e8b83c",
    },
    "source_repository": "https://github.com/ringsaturn/shachen/",
    "source_branch": "main",
    "source_directory": "docs/",
    "footer_icons": [
        {
            "name": "GitHub",
            "url": "https://github.com/ringsaturn/shachen",
            "html": (
                '<svg stroke="currentColor" fill="currentColor" stroke-width="0" '
                'viewBox="0 0 16 16"><path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 '
                "8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49"
                "-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01"
                "-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07"
                "-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12"
                "0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27s1.36.09 2 .27c1.53-1.04 2.2-.82"
                " 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65"
                " 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.012"
                ' 8.012 0 0 0 16 8c0-4.42-3.58-8-8-8z"></path></svg>'
            ),
            "class": "",
        },
    ],
}

# `0.2.1+g3fee821` in the sidebar says "ahead of the last tag" only to a reader
# who reads PEP 440 local versions for pleasure. Everyone else sees a version
# number and takes the page at its word. So say it in prose, on every page of
# the site built from main, with a way out to the version PyPI actually holds.
#
# The banner deliberately names no release number: between the version bump and
# the tag push, pyproject.toml already carries a version that has not shipped,
# and a banner that announced it would be pointing at a release that does not
# exist yet. `releases/latest` is always right. A build of a tagged commit gets
# no banner — it documents exactly the version it claims to.
#
# Keep the text to one short line: Furo lays the announcement out as a single
# row of the header's height, `white-space: nowrap`, so anything longer scrolls
# sideways instead of wrapping. This is HTML only — the PDF built from main
# carries the `+g<sha>` on its title page and nothing more.
if "+g" in release:
    html_theme_options["announcement"] = (
        f"You are reading the development docs, built from <code>main</code> "
        f"at <code>{release}</code> &mdash; ahead of the "
        f'<a href="https://github.com/ringsaturn/shachen/releases/latest">'
        f"latest release</a>."
    )

# -- LaTeX / PDF -------------------------------------------------------------

# xelatex, not the pdflatex default: the docs carry CJK (沙尘), the box-drawing
# characters of the pipeline diagram, and assorted arrows and math symbols,
# none of which pdflatex can typeset.
latex_engine = "xelatex"

latex_documents = [
    ("index", "shachen.tex", "shachen — DEBRA-Dust", author, "manual"),
]

# Sphinx indexes xelatex builds with xindy by default. MacTeX bundles it, the
# Debian texlive packages do not ship it at all, and installing it would drag
# clisp into CI for no gain: the index is a list of Python identifiers, ASCII
# to the last entry, which makeindex — part of texlive-binaries everywhere —
# sorts just as well. Turning it off keeps one indexer for every machine.
latex_use_xindy = False

latex_show_urls = "footnote"

# The cover image is not referenced from any page, so Sphinx would not copy it
# into the LaTeX build directory on its own.
latex_additional_files = ["img/social-preview.png"]

latex_elements = {
    "papersize": "a4paper",
    "pointsize": "11pt",
    # Both font families ship with TeX Live, so the PDF builds identically on
    # any machine — no system fonts to hunt down. (Debian and Ubuntu split the
    # font files into separate packages and attach them to the texlive ones as
    # Recommends, so a --no-install-recommends install has to name them; see
    # .github/actions/build-docs.)
    # DejaVu Sans Mono covers the box-drawing glyphs of the pipeline diagram,
    # Fandol covers the Chinese in the project name, and the newunicodechar
    # fallbacks catch the handful of symbols Latin Modern lacks in text mode
    # (verbatim is unaffected: there the mono font renders them directly).
    "preamble": r"""
\setmonofont{DejaVuSansMono.ttf}[
  BoldFont       = DejaVuSansMono-Bold.ttf,
  ItalicFont     = DejaVuSansMono-Oblique.ttf,
  BoldItalicFont = DejaVuSansMono-BoldOblique.ttf,
  Scale          = 0.85,
]
\usepackage{xeCJK}
\setCJKmainfont{FandolSong-Regular.otf}[BoldFont = FandolSong-Bold.otf]
% xeCJK claims the typographic dashes and quotes as CJK punctuation and pads
% them with CJK spacing, which mangles English prose ("the paper ' s"). The
% docs are English with two Chinese characters in them, so hand those code
% points back to the Latin font.
\xeCJKDeclareCharClass{Default}{"2013, "2014, "2018, "2019, "201C, "201D, "2026}
\usepackage{newunicodechar}
\newunicodechar{→}{\ensuremath{\rightarrow}}
\newunicodechar{≈}{\ensuremath{\approx}}
\newunicodechar{−}{\ensuremath{-}}
\newunicodechar{Δ}{\ensuremath{\Delta}}
\newunicodechar{µ}{\ensuremath{\mu}}
\newunicodechar{×}{\ensuremath{\times}}
""",
    # A cover page carrying the project banner ahead of the usual title page.
    # The banner is 2:1 and dark to the edges, so it goes full-bleed across the
    # paper width on a near-black page (the colour is sampled from the image's
    # own corners), which hides the seam above and below it.
    "maketitle": r"""
\begingroup
\thispagestyle{empty}
\definecolor{shachencover}{HTML}{080A0E}
\pagecolor{shachencover}
\null\vfill
\noindent\hspace*{-\dimexpr\oddsidemargin+1in\relax}%
\includegraphics[width=\paperwidth]{social-preview.png}%
\vfill\null
\newpage
\nopagecolor
\endgroup
\sphinxmaketitle
""",
}
