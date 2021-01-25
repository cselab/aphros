BUILD = build
MAKEFILE = $(BUILD)/Makefile

all: $(MAKEFILE)
	+make -C $(BUILD) $(target)
	+make -C $(BUILD) install

$(MAKEFILE):
	mkdir -p "$(BUILD)"
	(cd "$(BUILD)" && cmake ..)

install: $(MAKEFILE)
	+make -C $(BUILD) install

clean:
	rm -rf $(BUILD)

.PHONY: all clean install
