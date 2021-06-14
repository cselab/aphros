nx = 32
bsx = 32
m = $(nx) $(nx) 1
bs = $(bsx) $(bsx) 1
np = 1
tl = 600

include $(shell ap.makesim)

conf: b.dat

b.dat: gen
	./$<
