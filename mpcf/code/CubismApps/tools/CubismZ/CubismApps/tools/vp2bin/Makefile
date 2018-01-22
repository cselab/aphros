# module load hdf5-parallel/1.8.11
CC = mpic++
#CC = CC
bs = 32
hdf = 1
include ../../Makefile.config 

CPPFLAGS += -I../reader 
CPPFLAGS += $(extra)

#CSCS
#CPPFLAGS += -I${HDF5_INCLUDE_OPTS}

#ifeq "$(zlib)" "1"
#	MYFLAGS += -D_USE_ZLIB_
#endif
#
#ifeq "$(lz4)" "1"
#	MYFLAGS += -D_USE_LZ4_  
#endif    

all: vp2bin compbin testbin diffbin printbin

vp2bin: main.o WaveletCompressor.o
	$(CC) $(CPPFLAGS) WaveletCompressor.o  main.o  -o vp2bin $(LIBS) -ldl

main.o: main.cpp
	$(CC) $(CPPFLAGS) -c main.cpp -o main.o 
# -ldl
	
WaveletCompressor.o: ../../Compressor/source/WaveletCompressor.cpp
	$(CC) $(CPPFLAGS) -c ../../Compressor/source/WaveletCompressor.cpp -o WaveletCompressor.o

compbin: compbin.c
	$(CC) -DBS=$(bs) -O3 -o compbin compbin.c -lm

testbin: testbin.c
	$(CC) -DBS=$(bs) -O3 -o testbin testbin.c -lm	

diffbin: diffbin.c
	$(CC) -DBS=$(bs) -O3 -o diffbin diffbin.c -lm	

printbin: printbin.c
	$(CC) -DBS=$(bs) -O3 -o printbin printbin.c -lm	


clean:
	rm -f vp2bin compbin testbin diffbin printbin *.o

cleanall:
	rm -f vp2bin compbin testbin diffbin printbin *.o
