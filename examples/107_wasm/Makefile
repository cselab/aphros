BUILD = build
MAKEFILE = $(BUILD)/Makefile
CMAKE = emcmake cmake

all: $(MAKEFILE)
	+make -C $(BUILD) $(target)

configure: $(MAKEFILE)

$(MAKEFILE): Makefile CMakeLists.txt
	mkdir -p "$(BUILD)"
	cd "$(BUILD)" && $(CMAKE) .. -Daphrosjs_DIR="$${APHROS_PREFIX}/cmake"

clean:
	rm -rf $(BUILD)

aphros_doc.html: aphros_doc.md aphros_doc.css
	pandoc -s -c aphros_doc.css --metadata title='Aphros Explorer documentation' $< -o $@

run:
	cd "$(BUILD)" && emrun --serve_after_exit .

.PHONY: all clean configure run
