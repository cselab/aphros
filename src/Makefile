b = build

all:
	+make -C $(b) $(target)
	+make -C $(b) install 

install:
	+make -C $(b) install 

test:
	(cd build && ctest -j`nproc --all`)

.PHONY: all test install
