M = main
S = ../../../../src
T = libnorm_c.a

CH_CXXFLAGS = -I$S -std=c++11
$T: $M.o; ar r $@ $< && ranlib $@
.cpp.o:; $(CXX) -c $(CXXFLAGS) $(CH_CXXFLAGS) $< -o $@

.SUFFIXES: .cpp .o

.PHONY: clean
clean:; rm -f $M $M.o
