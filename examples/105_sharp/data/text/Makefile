default: blur.pgm
.PHONY: default

.tex.pdf:
	xelatex $<

.pdf.pgm:
	convert -density 300 $< -crop 256x256+0+0 -alpha remove $@

blur.pgm: main.pgm
	convert main.pgm -blur 8x2 blur.pgm

clean:
	$(RM) *.aux *.log *.pdf *.pgm


.SUFFIXES: .tex .pdf .pgm
