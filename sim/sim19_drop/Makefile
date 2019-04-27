
include $(shell ch.makesim)

clean::
	rm -f _m _bs _np _tl par.conf b.dat

mesh.conf: par.conf
mesh.conf: m ?= $(shell cat _m)
mesh.conf: bs ?= $(shell cat _bs)
mesh.conf: np ?= $(shell cat _np)
conf: tl ?= $(shell cat _tl)

par.conf: pargen
	./pargen
