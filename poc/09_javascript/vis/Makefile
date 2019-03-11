M4 = m4
M4FLAGS =

M = grid.html rect.html brush.html

all: $M
grid.html: grid.js
rect.html: rect.js

.in.html.html:
	@echo creating $@
	@trap 'rm -f $$t; exit 2' 1 2 3 15 && \
	t=/tmp/t.$$$$ && \
	$(M4) $(M4FLAGS) $< > $$t && \
	mv $$t $@

.SUFFIXES: .html .in.html

.PHONY: clean

clean:; rm -f $M
