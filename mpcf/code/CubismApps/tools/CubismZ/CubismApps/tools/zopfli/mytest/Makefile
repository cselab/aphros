all: mytest_zopfli mytest_deflate

mytest_zopfli: mytest.c
	gcc  -fopenmp -o mytest_zopfli mytest.c -I.. -L.. -lzopfli -lz -lm

mytest_deflate: mytest.c
	gcc -fopenmp -D_USE_DEFLATE_=1  -o mytest_deflate mytest.c -I.. -L.. -lzopfli -lz -lm

clean:
	rm -f mytest_zopfli mytest_deflate
