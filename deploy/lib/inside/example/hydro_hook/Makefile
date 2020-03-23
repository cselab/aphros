m = 32 32 32
bs = 16 16 16
np = 1
tl = 1440

B=build

default: runhook
.PHONY: default

$B/Makefile: CMakeLists.txt
	mkdir -p $B
	cd $B && ap.conf ..

$B/libhook.so: $B/Makefile hook.cpp
	make -C $B

runhook: conf
	LD_PRELOAD=build/libhook.so ap.run ap.mfer

add.conf:
	touch $@

conf: $B/libhook.so add.conf

include $(shell ap.makesim)

cleandat::
	rm -rf add.conf build
