all: d run

d: d.c div.h Makefile
	gcc -pedantic -Werror d.c -o d -lm

run: d
	./d

clean:
	rm -f d

.PHONY: run
