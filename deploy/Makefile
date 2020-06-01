DIRS = setenv cmake wrap eigen sim lib

all: $(DIRS)

$(DIRS):
	+make -C "$@"

.PHONY: $(DIRS) build
