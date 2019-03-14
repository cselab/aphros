S = ../../../..
include $S/make/js.mk

M = main curv

M4FLAGS = -I../../lib

all: $M

test: $M; atest test/*
.PHONY: clean test
clean:; rm -f $M
