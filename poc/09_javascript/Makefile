D = lib0/partstr lib0/matrix
T = vof_ends

install:
	for d in $D; \
	do (cd "$$d" && $(MAKE) $(MAKECMDGOALS)); \
	done

test:
	for d in $T; \
	do (cd "$$d" && $(MAKE) test); \
	done

.PHONY: install test
