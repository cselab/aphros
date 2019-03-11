S = ../../..
include $S/make/js.mk

M = main.js

all: $M
main.js: lib.js

run: $M; @nodejs main.js
.PHONY: clean run
clean:; rm -f $M
