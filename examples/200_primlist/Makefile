m = 64 64 64
bs = 32 32 32

default: sm_0000.vtk a.jpg

sm_0000.vtk: b.dat Makefile std.conf
	make cleanrun

a.jpg: sm_0000.vtk
	which pvba1tch && pvbatch ./plot.py || touch a.jpg

.PHONY: default

include $(shell ap.makesim)

clean::
	rm -f a.jpg
