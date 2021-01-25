default: blur.png
.PHONY: default

.tex.pdf:
	xelatex $<

.pdf.png:
	convert -density 300 $< -crop 256x256+0+0 -alpha remove $@

blur.png: main.png
	convert main.png -blur 8x2 blur.png

clean:
	$(RM) *.aux *.log *.pdf *.png


.SUFFIXES: .tex .pdf .png
