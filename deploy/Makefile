BUILD = build
MAKEFILE = $(BUILD)/Makefile
CMAKE = cmake

all: $(MAKEFILE)
	+make -C $(BUILD) $(target)
	+make -C $(BUILD) install

$(MAKEFILE):
	mkdir -p "$(BUILD)"
	cd "$(BUILD)" && $(CMAKE) ..

install: $(MAKEFILE)
	+make -C $(BUILD) install

clean:
	rm -rf $(BUILD)

.PHONY: all clean install
