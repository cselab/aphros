M = main

$M: $M.o Makefile

CH_CXXFLAGS = -I$(CHPREFIX)/include/hydro
CH_LDFLAGS = -L$(CHPREFIX)/lib -ldistr -lnormal -lsuspender

$M: $M.o; $(CXX) $< $(CH_LDFLAGS) -o $@
.cpp.o:; $(CXX) -c $(CXXFLAGS) $(CH_CXXFLAGS) $< -o $@

.SUFFIXES: .cpp .o

.PHONY: clean
clean:; rm -f $M $M.o
