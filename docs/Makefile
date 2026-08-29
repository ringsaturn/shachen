# Minimal makefile for Sphinx documentation.
# `make html` then open _build/html/index.html.

SPHINXOPTS    ?= -W --keep-going
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# PDF needs a TeX distribution (MacTeX on macOS, texlive-xetex + texlive-latex-extra
# + texlive-fonts-extra + latexmk on Debian/Ubuntu). Everything it typesets with
# ships with TeX Live, so no fonts to install:
#     make -C docs latexpdf   ->  docs/_build/latex/shachen.pdf
