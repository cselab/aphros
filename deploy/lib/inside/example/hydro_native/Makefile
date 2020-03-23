m = 32 32 32
bs = 16 16 16
np = 1
tl = 1440

default: cleanrun
.PHONY: default

add.conf:
	touch $@

conf: add.conf

include $(shell ap.makesim)

cleandat::
	rm -rf add.conf
