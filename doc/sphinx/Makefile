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

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile all install web default man1 man3

clean:
	rm -rf build

%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
