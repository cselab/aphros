D = lib0/m4 lib0/partstr lib0/matrix lib0/vof
T = example/vof_ends example/cell_ends example/full

install:
	for d in $D; \
	do (cd "$$d" && $(MAKE) install); \
	done

test:
	for d in $T; \
	do (cd "$$d" && $(MAKE) test); \
	done

.PHONY: install test
