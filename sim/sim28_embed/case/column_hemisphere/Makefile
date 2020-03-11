m = 32 64 32
np = 8
#m = 64 128 64
#m = 128 256 128
#bs = 32 32 32
tl = 600

default:
	killall ap.mfer || true
	make cleanall
	nice make run

add.conf: add.conf.gen
	./add.conf.gen

conf: add.conf

include $(shell ap.makesim)

clean::
	rm -vrf sc
