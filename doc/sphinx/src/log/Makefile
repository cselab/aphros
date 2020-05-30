PANDOC = pandoc
T = log

all: $T.html

$T.html: $T.rst Makefile
	$(PANDOC) -s -c main.css --mathjax -o $@ $<
