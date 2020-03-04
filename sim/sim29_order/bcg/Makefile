m = 128 128 1
bs = 16 16 1
np = 8
tl = 1440

B=build

default: conf
.PHONY: default

$B/Makefile: CMakeLists.txt
	mkdir -p $B
	cd $B && ap.conf ..

$B/libhook.so: $B/Makefile hook.cpp
	make -C $B

conf: $B/libhook.so

include $(shell ap.makesim)

clean::
	rm -f omz_*.dat par.conf
