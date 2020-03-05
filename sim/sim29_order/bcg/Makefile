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

par.conf:
	nx=`cut -d' ' -f1 <<< "$m"` ; \
	dt=`bc -l <<< "0.5 / $$nx"` ; \
	echo set double dt0 $$dt > par.conf
.PHONY: par.conf

runhook: conf
	LD_PRELOAD=build/libhook.so ap.run ap.mfer

conf: $B/libhook.so par.conf

include $(shell ap.makesim)

clean::
	rm -f omz_*.dat par.conf
