CFLAGS=-O2 -g
INCLUDE_PATH=../include 

all:
	g++ -I${INCLUDE_PATH} ${CFLAGS} -c wavelets.c haar.c
	mv *o ../lib

test:
	g++ -I${INCLUDE_PATH} ${CFLAGS} -DTEST=1 wavelets.c haar.c -o wavelets ../src/utils.c ../src/bitstream.c -lgsl -lgslcblas

clean:
	rm -f *o
