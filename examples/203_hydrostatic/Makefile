m = 128 128 1
np = 1
bs = 32 32 1

default: help
.PHONY: default

include $(shell ap.makesim)

single:
	> add.conf
	make m="128 128 1" cleanrun

multi:
	cp multi.conf add.conf
	make m="32 32 1" cleanrun

.PHONY: single multi
