m = 32 32 1
#m = 64 64 1
#m = 128 32 1
np = 1
#m = 512 128 1
#np = 256
#m = 256 64 1
#np = 64
bs = 32 32 1
tl = 1440

default:
	killall ap-mfer || true
	make cleanall
	make run

include $(shell ap-makesim)

clean::
	rm -vrf sc
