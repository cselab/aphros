CC = gcc
CFLAGS = -std=c99 -Wall -pedantic -g
CXX = g++
CXXFLAGS =

# uncomment to disable array headers
#CFLAGS += -DWITHOUT_HEADER

testfpzip: testfpzip.o ../lib/libfpzip.a
	$(CXX) $(CXXFLAGS) -L../lib testfpzip.o -lfpzip -o $@

testfpzip.o: testfpzip.c
	$(CC) $(CFLAGS) -I../inc -c $<

clean:
	@rm -f testfpzip testfpzip.o
