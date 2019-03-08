M = main
S = ../../../src

$M: $M.o

CH_CXXFLAGS = -I $S -std=c++11

$M: $M.o; $(CXX) $< $(LDFLAGS) -o $@
.cpp.o:; $(CXX) -c $(CXXFLAGS) $(CH_CXXFLAGS) $< -o $@

.SUFFIXES: .cpp .o

.PHONY: clean
clean:; rm -f $M $M.o
