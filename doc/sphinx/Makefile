SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SOURCEDIR     = src
BUILDDIR      = build

default: install

all: install web

web: singlehtml

man1 man3 man7: man
	s=$(subst man,,$@) ; \
	d="$(CHPREFIX)/man/man$$s" ; \
	mkdir -p "$$d" ; \
	cp -v build/man/*.$$s "$$d/"

install: man1 man3 man7
	mandb -q --user-db

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile all install web default man1 man3 man7

clean:
	rm -rf build

%: Makefile
	$(SPHINXBUILD) $(SPHINXOPTS) -b $@ "$(SOURCEDIR)" "$(BUILDDIR)/$@"
