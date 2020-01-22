m = 128 64 1
bs = 8 8 1
np = 128
tl = 1440

default:
	killall ap-mfer || true
	make cleanall
	make run

include $(shell ap-makesim)

clean::
	rm -vrf sc
