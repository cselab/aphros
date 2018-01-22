CC = gcc
CFLAGS = -O3

all: libfpc.a main

libfpc.a: fpc.o
	ar rs libfpc.a fpc.o
	
fpc.o: fpc.c
	$(CC) $(CFLAGS) -c fpc.c 

main: main.c
	$(CC) $(CFLAGS) -o main main.c -L. -lfpc

clean:
	rm -f libfpc.a fpc.o main.o main
