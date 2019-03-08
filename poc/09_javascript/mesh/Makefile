M = main
S = ../../../src

$M: $M.o

CH_CXXFLAGS = -I$S -std=c++11
CH_LDFLAGS = -L$S/build/solver -L$S/build/util -lsolver -lnormal -lsuspender

$M: $M.o; $(CXX) $< $(CH_LDFLAGS) -o $@
.cpp.o:; $(CXX) -c $(CXXFLAGS) $(CH_CXXFLAGS) $< -o $@

.SUFFIXES: .cpp .o

.PHONY: clean
clean:; rm -f $M $M.o
