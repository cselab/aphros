BUILD = build
MAKEFILE = $(BUILD)/Makefile
CMAKE = cmake

all: $(MAKEFILE)
	+make -C $(BUILD) $(target)
	+make -C $(BUILD) install 

$(MAKEFILE):
	mkdir -p "$(BUILD)"
	(cd "$(BUILD)" && $(CMAKE) ..)

install: $(MAKEFILE)
	+make -C $(BUILD) install 

test:
	(cd $(BUILD) && ctest -j`nproc --all`)

.PHONY: all install test
