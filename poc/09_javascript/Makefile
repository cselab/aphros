M = main

$M: $M.o

CH_CXXFLAGS = -I../../src -std=c++11
CH_LDFLAGS = -L../../src/build/util -lsuspender

$M: $M.o; $(CXX) $< $(CH_LDFLAGS) -o $@
.cpp.o:; $(CXX) -c $(CXXFLAGS) $(CH_CXXFLAGS) $< -o $@

.SUFFIXES: .cpp .o

.PHONY: clean
clean:; rm -f $M $M.o
