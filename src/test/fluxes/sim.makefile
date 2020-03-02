m = 16 16 1
bs = 16 16 1
np = 1
tl = 10

include $(shell ap.makesim)

clean::
	rm -f vx add.conf refpath
