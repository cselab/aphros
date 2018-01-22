CFLAGS=-O2 -g

all:
	g++ ${CFLAGS} -o file_compress -DCOMPRESS=1 file_compress.c -I../../include -L../../lib -lgsl -lz -lgslcblas -lisabela
	g++ ${CFLAGS} -o file_decompress file_compress.c -I../../include -L../../lib -lgsl -lz -lgslcblas -lisabela

clean:
	rm -f *o;
	rm -f file_compress file_decompress;
