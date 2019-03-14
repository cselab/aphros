S = ../../..
PREFIX = ..
include $S/make/html.mk

M = grid.html rect.html brush.html
M4FLAGS = -I$(PREFIX)/lib

all: $M
grid.html: grid.js
rect.html: rect.js
.PHONY: clean
clean:; rm -f $M
