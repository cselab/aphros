.SUFFIXES:
SRC = ../../src
WRK = .work
CXX = emcc
CXXFLAGS = -Wall -O3 -DNDEBUG
WASMFLAGS = -s ALLOW_MEMORY_GROWTH=1
ALL = hydro.js electrochem.js aphros.js
1 = -D_USE_MPI_=0
all: $(ALL)
# $(WRK)/libaphros_js.a:; $(MAKE) -f $(SRC)/Makefile_legacy 'WRK=$(WRK)' 'SRC=$(SRC)' 'CXX=$(CXX)' 'CC=$(CXX)' 'CXXFLAGS=$(CXXFLAGS)' $@
L = ../../src/build/libaphros.a
hydro.js: hydro.cpp 
	$(CXX) $< $(CXXFLAGS) $(WASMFLAGS) -I$(SRC) $L $1 \
		-s EXPORTED_FUNCTIONS='["_main", "_AddVelocityAngle", "_SetRuntimeConfig", "_SetExtraConfig", "_TogglePause", "_Init", "_Spawn", "_SetMesh" ,"_GetLines", "_SetCoal"]' \
		-s EXPORTED_RUNTIME_METHODS='["cwrap", "ccall"]' -o $@
electrochem.js: electrochem.cpp $L
	$(CXX) $< $(CXXFLAGS) $(WASMFLAGS) -I$(SRC) $L $1 \
		-s EXPORTED_FUNCTIONS='["_main", "_SetRuntimeConfig", "_SetExtraConfig", "_TogglePause", "_Init", "_Spawn", "_SetMesh" ,"_GetLines"]' \
		-s EXPORTED_RUNTIME_METHODS='["cwrap", "ccall"]' -o $@

aphros.js: aphros.cpp $L conf/base.conf conf/std.conf
	$(CXX) $< $(CXXFLAGS) $(WASMFLAGS) -I$(SRC) $L $1 \
		--preload-file conf \
		-s EXPORTED_FUNCTIONS='["_main", "_GetConfigString", "_GetConfigDouble", "_SetRuntimeConfig", "_SetExtraConfig", "_TogglePause", "_Spawn", "_SetMesh" ,"_GetLines"]' \
		-s EXPORTED_RUNTIME_METHODS='["cwrap", "ccall"]' -o $@

conf/base.conf: ../../deploy/scripts/sim_base.conf
	cp $< $@

aphros_doc.html: aphros_doc.md aphros_doc.css
	pandoc -s -c aphros_doc.css --metadata title='Aphros Explorer documentation' $< -o $@

diffusion.js: diffusion.cpp $L
	$(CXX) $< $(CXXFLAGS) $(WASMFLAGS) -I$(SRC) $L $1 \
		-s EXPORTED_FUNCTIONS='["_main", "_MulDiffusion", "_AddVelocityAngle", "_TogglePause", "_Init", "_SetMesh"]' \
		-s EXPORTED_RUNTIME_METHODS='["cwrap", "ccall"]' -o $@
clean:; rm -vf $(ALL) conf/base.conf && $(MAKE) -f $(SRC)/Makefile_legacy 'WRK=$(WRK)' 'SRC=$(SRC)' clean
.PHONY: clean all
