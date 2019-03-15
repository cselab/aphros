S = ../../../..
include $S/make/js.mk

M = main
M4FLAGS = -I../../lib

all: $M
lint: $M
	jshint --reporter=unix --verbose -c data/conf.json $M

test: $M; atest test/*
.PHONY: clean test lint
clean:; rm -f $M
