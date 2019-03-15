S = ../../..
PREFIX = ..
include $S/make/html.mk
include $S/make/js.mk

M = grid.html rect.html brush.html g.js
M4FLAGS = -I$(PREFIX)/lib

all: $M
grid.html: grid.js
rect.html: rect.js

lint: g.js
	jshint --reporter=unix --verbose -c data/conf.json $^

.PHONY: clean all lint

clean:; rm -f $M
