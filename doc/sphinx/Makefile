.SUFFIXES:
.PHONY: all clean default install latexpdf man1 man3 man7 singlehtml singlehtml_data web
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

install: man1 man3 man7; mandb -q --user-db
clean:; rm -rf build
include img_pdf.mk
latexpdf: $(IMG_PDF)
%.pdf: %.svg; inkscape --export-filename=$@ $< || inkscape --export-pdf=$@ $<
latexpdf: Makefile; $(SPHINXBUILD) $(SPHINXOPTS) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)/$@"
man: Makefile; $(SPHINXBUILD) $(SPHINXOPTS) -b $@ "$(SOURCEDIR)" "$(BUILDDIR)/$@"
singlehtml: Makefile; $(SPHINXBUILD) $(SPHINXOPTS) -b $@ "$(SOURCEDIR)" "$(BUILDDIR)/$@"
