m = 128 128 128
np = 128
bs = 16 16 16
tl = 600

default:
	killall ap.mfer || true
	make cleanall
	make run

include $(shell ap.makesim)

clean::
	rm -vrf sc
