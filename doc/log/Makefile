PANDOC = pandoc
T = log

all: $T.html

$T.html: $T.md Makefile
	$(PANDOC) -s -c main.css --mathjax -o $@ $<
