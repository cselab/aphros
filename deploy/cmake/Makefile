all:
	mkdir -p build
	( cd build && . ap.setenv && cmake ../.. -DCMAKE_INSTALL_PREFIX=$$CHPREFIX )
	( cd build && make && make install )

.PHONY: all
