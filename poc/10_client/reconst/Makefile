M = main

$M: $M.o

CH_CXXFLAGS = -I$(CHPREFIX)/include/hydro

$M: $M.o; $(CXX) $< $(LDFLAGS) -o $@
.cpp.o:; $(CXX) -c $(CXXFLAGS) $(CH_CXXFLAGS) $< -o $@

.SUFFIXES: .cpp .o

.PHONY: clean
clean:; rm -f $M $M.o
