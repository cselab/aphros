.SUFFIXES:
SRC = ../../src
WRK = .work
CXX = emcc
CXXFLAGS = -Wall -O3
WASMFLAGS = -s ALLOW_MEMORY_GROWTH=1
APHROS_PREFIX ?= $(HOME)/hy
ALL = hydro.js
1 = -D_USE_MPI_=0
all: $(ALL)
$(WRK)/libaphros_js.a:; $(MAKE) -f $(SRC)/Makefile_legacy 'WRK = $(WRK)' 'SRC = $(SRC)' 'CXX = $(CXX)' 'CC = $(CXX)' 'CXXFLAGS = $(CXXFLAGS)' $@
hydro.js: hydro.cpp $(WRK)/libaphros_js.a
	$(CXX) $< $(CXXFLAGS) $(WASMFLAGS) -I$(SRC) $(WRK)/libaphros_js.a $1 \
		-s EXPORTED_FUNCTIONS='["_main", "_AddVelocityAngle", "_SetExtraConfig", "_TogglePause", "_Init", "_Spawn", "_SetMesh" ,"_GetLines", "_SetCoal"]' \
		-s EXPORTED_RUNTIME_METHODS='["cwrap", "ccall"]' -o $@
clean:; rm -vf $(ALL)
.PHONY: clean all
