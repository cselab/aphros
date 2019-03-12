S = ../../../..
include $S/make/js.mk

M = main.js
NODE = node

M4FLAGS = -I../lib

all: $M

run: $M; @$(NODE) $M 0.1 1 0.3
.PHONY: clean run
clean:; rm -f $M
