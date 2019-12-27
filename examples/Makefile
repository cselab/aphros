DIRS = $(wildcard ???_*)
OUT = /dev/null

all: $(DIRS)

$(DIRS):
	@echo "$@"
	@(cd "$@" && ap-conf .) > $(OUT)
	@+make -s -C "$@" > $(OUT)

clean:
	for d in $(DIRS) ; do make -C "$$d" clean ; done

.PHONY: all clean $(DIRS)
