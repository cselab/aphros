all: mytest_zstd 

mytest_zstd: mytest.c
	gcc  -fopenmp -o mytest_zstd mytest.c -I../lib  -L.  -lzstd  -lm

clean:
	rm -f mytest_zstd
