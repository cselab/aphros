SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SOURCEDIR     = src
BUILDDIR      = build

default: install

all: install web

web: singlehtml singlehtml_data

man1 man3 man7: man
	s=$(subst man,,$@) ; \
	d="$(APHROS_PREFIX)/man/man$$s" ; \
	mkdir -p "$$d" ; \
	cp -v build/man/*.$$s "$$d/"

singlehtml_data:
	cp -r src/log/log??_* build/singlehtml/
	rsync -iaR .././dashboard/*.{js,html,css} build/singlehtml/

install: man1 man3 man7
	mandb -q --user-db

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile all install web default man1 man3 man7 singlehtml_data

clean:
	rm -rf build

%: Makefile
	$(SPHINXBUILD) $(SPHINXOPTS) -b $@ "$(SOURCEDIR)" "$(BUILDDIR)/$@"
