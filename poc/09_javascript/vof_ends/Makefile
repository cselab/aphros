S = ../../..
include $S/make/js.mk

M = main.js
NODE = node

M4FLAGS = -I../lib

all: $M

run: $M; @$(NODE) $M data/ij
test: $M; atest test/*
.PHONY: clean run test
clean:; rm -f $M
