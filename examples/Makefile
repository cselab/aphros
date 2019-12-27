DIRS = $(wildcard ???_*)
OUT = /dev/null

all:
	+make build
	+make run

build: $(DIRS)

$(DIRS):
	@echo "$@"
	@(cd "$@" && ap-conf .) > $(OUT)
	@+make -s -C "$@" > $(OUT)

run:
	@for d in $(DIRS) ; do echo "$$d" ; (cd "$$d" && ./run) ; done

clean:
	@for d in $(DIRS) ; do echo "$$d" ; make -C "$$d" clean ; done

.PHONY: all build run clean $(DIRS)
