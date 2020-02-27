#m = 64 128 64
m = 128 256 128
np = 128
bs = 16 16 16
tl = 600

default:
	killall ap.mfer || true
	make cleanall
	make run

add.conf: add.conf.gen
	./add.conf.gen

conf: add.conf

include $(shell ap.makesim)

clean::
	rm -vrf sc
