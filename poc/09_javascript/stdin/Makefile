S = ../../..
include $S/make/js.mk

M = main.js
NODE = node

all: $M

run: $M; @$(NODE) $M
.PHONY: clean run
clean:; rm -f $M
