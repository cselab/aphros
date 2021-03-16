.SUFFIXES:
SRC = ../../src
WRK = .work
CXX = emcc
CXXFLAGS = -Wall -O3 -DNDEBUG
WASMFLAGS = -s ALLOW_MEMORY_GROWTH=1
ALL = hydro.js electrochem.js aphros.js
1 = -D_USE_MPI_=0
all: $(ALL)
$(WRK)/libaphros_js.a:; $(MAKE) -f $(SRC)/Makefile_legacy 'WRK=$(WRK)' 'SRC=$(SRC)' 'CXX=$(CXX)' 'CC=$(CXX)' 'CXXFLAGS=$(CXXFLAGS)' $@
hydro.js: hydro.cpp $(WRK)/libaphros_js.a
	$(CXX) $< $(CXXFLAGS) $(WASMFLAGS) -I$(SRC) $(WRK)/libaphros_js.a $1 \
		-s EXPORTED_FUNCTIONS='["_main", "_AddVelocityAngle", "_SetRuntimeConfig", "_SetExtraConfig", "_TogglePause", "_Init", "_Spawn", "_SetMesh" ,"_GetLines", "_SetCoal"]' \
		-s EXPORTED_RUNTIME_METHODS='["cwrap", "ccall"]' -o $@
electrochem.js: electrochem.cpp $(WRK)/libaphros_js.a
	$(CXX) $< $(CXXFLAGS) $(WASMFLAGS) -I$(SRC) $(WRK)/libaphros_js.a $1 \
		-s EXPORTED_FUNCTIONS='["_main", "_SetRuntimeConfig", "_SetExtraConfig", "_TogglePause", "_Init", "_Spawn", "_SetMesh" ,"_GetLines"]' \
		-s EXPORTED_RUNTIME_METHODS='["cwrap", "ccall"]' -o $@

aphros.js: aphros.cpp $(WRK)/libaphros_js.a
	$(CXX) $< $(CXXFLAGS) $(WASMFLAGS) -I$(SRC) $(WRK)/libaphros_js.a $1 \
		--preload-file conf \
		-s EXPORTED_FUNCTIONS='["_main", "_GetConfigString", "_GetConfigDouble", "_SetRuntimeConfig", "_SetExtraConfig", "_TogglePause", "_Spawn", "_SetMesh" ,"_GetLines"]' \
		-s EXPORTED_RUNTIME_METHODS='["cwrap", "ccall"]' -o $@

aphros_doc.html: aphros.md
	pandoc -s -c pandoc.css $< -o $@

diffusion.js: diffusion.cpp $(WRK)/libaphros_js.a
	$(CXX) $< $(CXXFLAGS) $(WASMFLAGS) -I$(SRC) $(WRK)/libaphros_js.a $1 \
		-s EXPORTED_FUNCTIONS='["_main", "_MulDiffusion", "_AddVelocityAngle", "_TogglePause", "_Init", "_SetMesh"]' \
		-s EXPORTED_RUNTIME_METHODS='["cwrap", "ccall"]' -o $@
clean:; rm -vf $(ALL) && $(MAKE) -f $(SRC)/Makefile_legacy 'WRK=$(WRK)' 'SRC=$(SRC)' clean
.PHONY: clean all
