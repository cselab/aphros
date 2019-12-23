# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = python -msphinx
SPHINXPROJ    = hydro
SOURCEDIR     = src
BUILDDIR      = build

default: install

all: install web

web: singlehtml

man1: man
	s=1 ; \
	d="$(CHPREFIX)/man/man$$s" ; \
	mkdir -p "$$d" ; \
	cp -v build/man/*.$$s "$$d/"

man3: man
	s=3 ; \
	d="$(CHPREFIX)/man/man$$s" ; \
	mkdir -p "$$d" ; \
	cp -v build/man/*.$$s "$$d/"

install: man1 man3
	mandb --user-db

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile all install web default man1 man3

clean:
	rm -rf build

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
