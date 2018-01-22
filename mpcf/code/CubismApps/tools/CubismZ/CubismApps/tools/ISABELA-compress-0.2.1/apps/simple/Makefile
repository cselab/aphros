CFLAGS=-O2 -g

all:
	g++ ${CFLAGS} -o simple simple.c -I../../include -L../../lib -lgsl -lz -lgslcblas -lisabela

clean:
	rm -f *o;
	rm -f file_compress file_decompress;
