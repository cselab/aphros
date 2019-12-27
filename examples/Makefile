DIRS = $(wildcard ???_*)
VERBOSE =
OUT = $(if $(VERBOSE),/dev/stdout,/dev/null)

all:
	+make build
	+make run

build: $(DIRS)

$(DIRS):
	@echo "$@"
	@(cd "$@" && ap-conf .) >$(OUT) 2>&1
	@+make -s -C "$@" >$(OUT) 2>&1

run:
	@for d in $(DIRS) ; do echo "$$d" ; (cd "$$d" && ./run >$(OUT) 2>&1 ) ; done

clean:
	@for d in $(DIRS) ; do echo "$$d" ; make -C "$$d" clean >$(OUT) 2>&1 ; done

.PHONY: all build run clean $(DIRS)
