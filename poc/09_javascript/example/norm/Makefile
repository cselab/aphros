S = ../../../..
include $S/make/js.mk

M = main.js
NODE = node

M4FLAGS = -I../lib

all: $M
main.js: lib.js

run: $M; @$(NODE) $M
.PHONY: clean run
clean:; rm -f $M
