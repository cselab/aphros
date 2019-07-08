b = build

all:
	+make -C $(b)
	+make -C $(b) install 

install:
	+make -C $(b) install 

test:
	(cd build && ctest -j`nproc`)

.PHONY: all test install
