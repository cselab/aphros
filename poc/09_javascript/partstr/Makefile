S = ../../..
PREFIX = ..

include $S/make/js.mk

M = partstr.js
NODE = node

all: $M

run: $M; @$(NODE) $M
install: $M
	mkdir -p $(DESTDIR)$(PREFIX)/lib
	cp $M $(PREFIX)/lib

.PHONY: clean run install
clean:; rm -f $M
