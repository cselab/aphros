m = 192 96 96
bs = 16 16 16
np = 72
tl = 10

include $(shell ch.makesim)

clean::
	rm -f _m _bs _np _tl

run1:
	make m="64 32 32" np=1 run

run72:
	make m="192 96 96" np=72 run

submit1:
	make m="64 32 32" np=1 submit

submit72:
	make m="192 96 96" np=72 submit


.PHONY: run1 run72
