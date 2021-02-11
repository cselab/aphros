.SUFFIXES:
S = ../../src
B = .work
CXX = emcc
CXXFLAGS = -Wall -O3
WASMFLAGS = -s ALLOW_MEMORY_GROWTH=1
APHROS_PREFIX ?= $(HOME)/hy
ALL = hydro.js
1 = -D_USE_MPI_=0
all: $(ALL)
$B/libaphros_js.a:; $(MAKE) -f $S/Makefile_legacy 'B = $B' 'S = $S' 'CXX = $(CXX)' 'CC = $(CXX)' 'CXXFLAGS = $(CXXFLAGS)' $@
hydro.js: hydro.cpp $B/libaphros_js.a
	$(CXX) $< $(CXXFLAGS) $(WASMFLAGS) -I$S $B/laphros_js.a $1 \
		-s EXPORTED_FUNCTIONS='["_main", "_AddVelocityAngle", "_SetExtraConfig", "_TogglePause", "_Init", "_SetMesh"]' \
		-s EXPORTED_RUNTIME_METHODS='["cwrap", "ccall"]' -o $@
clean:; rm -vf $(ALL)
.PHONY: clean all
