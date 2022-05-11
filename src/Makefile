BUILD = build
MAKEFILE = $(BUILD)/Makefile
CMAKE = cmake

JS_BUILD = buildjs
JS_MAKEFILE = $(JS_BUILD)/Makefile
JS_CMAKE = emcmake cmake
JS_FLAGS = -DUSE_OPENMP=0 -DUSE_TESTS=0 -DUSE_BACKEND_CUBISM=0 -DUSE_BACKEND_LOCAL=0 -DUSE_EXPLORER=0 -DUSE_AVX=0 -DUSE_HYPRE=0 -DUSE_MPI=0 -DUSE_HDF=0 -DUSE_JS=1 -DUSE_MFER=0 -DUSE_CONF2PY=0

all: $(MAKEFILE)
	+make -C $(BUILD) $(target)
	+make -C $(BUILD) install

$(MAKEFILE): Makefile
	mkdir -p "$(BUILD)"
	cd "$(BUILD)" && $(CMAKE) ..

install: $(MAKEFILE)
	+make -C $(BUILD) install

$(JS_MAKEFILE): Makefile
	mkdir -p "$(JS_BUILD)"
	cd "$(JS_BUILD)" && $(JS_CMAKE) .. $(JS_FLAGS)

js: $(JS_MAKEFILE)
	+make -C $(JS_BUILD) $(target)
	+make -C $(JS_BUILD) install

test:
	cd $(BUILD) && ctest -j`nproc --all`

clean:
	rm -rf $(BUILD)

.PHONY: all clean install test
