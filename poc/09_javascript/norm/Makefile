S = ../../..
include $S/make/js.mk

M = main.js
NODE = node

all: $M
main.js: lib.js

run: $M; @$(NODE) $M.js
.PHONY: clean run
clean:; rm -f $M
