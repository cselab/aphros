bs ?= 16

MYFLAGS =  -D_BLOCKSIZE_=$(bs) -D_FLOAT_PRECISION_ -D_SP_COMP_ -g -fopenmp -I../../../Cubism/source/

mytest: main.o WaveletCompressor.o
	$(CC) WaveletCompressor.o  main.o  -o mytest  $(MYFLAGS) -lz -lstdc++

main.o: main.cpp
	$(CC) $(MYFLAGS) -c main.cpp -o main.o
	
WaveletCompressor.o: ../../MPCFnode/source/WaveletCompressor.cpp
	$(CC) $(MYFLAGS) -c ../../MPCFnode/source/WaveletCompressor.cpp -o WaveletCompressor.o
	
clean: *.o
	rm mytest *.o