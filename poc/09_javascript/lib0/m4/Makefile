PREFIX = ../..

M = mh.m4
NODE = node

install: $M
	mkdir -p $(DESTDIR)$(PREFIX)/lib && \
	cp $< $(PREFIX)/lib

.PHONY: clean run install
clean:; rm -f $M
