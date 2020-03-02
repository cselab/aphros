m = 8 16 1
bs = 8 8 1
np = 1
tl = 10

include $(shell ap.makesim)

clean::
	rm -f vx add.conf refpath
