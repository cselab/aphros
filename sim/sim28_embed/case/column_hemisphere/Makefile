#m = 64 128 64
m = 128 256 128
np = 128
bs = 16 16 16
tl = 600

default:
	./add.conf.gen
	killall ap.mfer || true
	make cleanall
	make run

include $(shell ap.makesim)

clean::
	rm -vrf sc
