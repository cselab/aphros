S = ../../..
include $S/make/html.mk

M = grid.html rect.html brush.html
M4FLAGS = -I../lib0

all: $M
grid.html: grid.js
rect.html: rect.js
.PHONY: clean
clean:; rm -f $M
