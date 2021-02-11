.SUFFIXES:
CXX = emcc
CXXFLAGS = -Wall -O3
WASMFLAGS = -s ALLOW_MEMORY_GROWTH=1
APHROS_PREFIX ?= $(HOME)/hy
ALL = hydro.js
1 = -D_USE_MPI_=0
all: $(ALL)
hydro.js: hydro.cpp
	$(CXX) $< $(CXXFLAGS) $(WASMFLAGS) -I$(APHROS_PREFIX)/include/aphros -L$(APHROS_PREFIX)/lib -laphros_js $1 \
		-s EXPORTED_FUNCTIONS='["_main", "_AddVelocityAngle", "_SetExtraConfig", "_TogglePause", "_Init", "_SetMesh"]' \
		-s EXPORTED_RUNTIME_METHODS='["cwrap", "ccall"]' -o $@
clean:; rm -vf $(ALL)
.PHONY: clean all
