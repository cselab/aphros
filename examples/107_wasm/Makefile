.SUFFIXES:
S = ../../src
CXX = emcc
CXXFLAGS = -Wall -O3
WASMFLAGS = -s ALLOW_MEMORY_GROWTH=1
APHROS_PREFIX ?= $(HOME)/hy
ALL = hydro.js
1 = -D_USE_MPI_=0
all: $(ALL)
libaphros_js.a:; $(MAKE) -f ../../src/Makefile_legacy 'S = $(S)' 'CXX = $(CXX)' 'CC = $(CXX)' 'CXXFLAGS = $(CXXFLAGS)' libaphros_js.a
hydro.js: hydro.cpp libaphros_js.a
	$(CXX) $< $(CXXFLAGS) $(WASMFLAGS) -I$S laphros_js.a $1 \
		-s EXPORTED_FUNCTIONS='["_main", "_AddVelocityAngle", "_SetExtraConfig", "_TogglePause", "_Init", "_SetMesh"]' \
		-s EXPORTED_RUNTIME_METHODS='["cwrap", "ccall"]' -o $@
clean:; rm -vf $(ALL)
.PHONY: clean all
