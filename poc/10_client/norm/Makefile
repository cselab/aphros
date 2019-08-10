M = main

$M: $M.o

CH_CXXFLAGS = -I$(CHPREFIX)/include/hydro
CH_LDFLAGS = -L$(CHPREFIX)/lib -lhydro

$M: $M.o; $(CXX) $< $(CH_LDFLAGS) -o $@
.cpp.o:; $(CXX) -c $(CXXFLAGS) $(CH_CXXFLAGS) $< -o $@

.SUFFIXES: .cpp .o

.PHONY: clean
clean:; rm -f $M $M.o
