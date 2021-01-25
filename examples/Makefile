DIRS = $(wildcard ???_*)
VERBOSE =
REDIR = $(if $(VERBOSE),,>/dev/null 2>&1)

all:
	+make build
	+make run

build: $(DIRS)

$(DIRS):
	@echo "$@"
	@(cd "$@" && ! test -f CMakeLists.txt ||  ap.conf . ) $(REDIR)
	@+make -s -C "$@" $(REDIR)

run:
	@for d in $(DIRS) ; do echo "$$d" ; (cd "$$d" && ./run $(REDIR) ) ; done

clean:
	@for d in $(DIRS) ; do echo "$$d" ; make -C "$$d" clean $(REDIR) ; done

.PHONY: all build run clean $(DIRS)
