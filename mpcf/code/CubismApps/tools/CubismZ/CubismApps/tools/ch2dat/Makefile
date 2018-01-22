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

all: ch2dat compdat

ch2dat: main.o WaveletCompressor.o
	$(CC) $(CPPFLAGS) WaveletCompressor.o  main.o  -o ch2dat $(LIBS) -ldl

main.o: main.cpp
	$(CC) $(CPPFLAGS) -c main.cpp -o main.o 
# -ldl


WaveletCompressor.o: ../../Compressor/source/WaveletCompressor.cpp
	$(CC) $(CPPFLAGS) -c ../../Compressor/source/WaveletCompressor.cpp -o WaveletCompressor.o

compdat: compdat.c
	$(CC) $(CPPFLAGS) compdat.c -o compdat
	
clean:
	rm -f ch2dat compdat *.o

cleanall:
	rm -f ch2dat compdat *.o
