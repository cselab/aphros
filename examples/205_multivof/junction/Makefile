nx = 64
bsx = 64
m = $(nx) $(nx) 1
bs = $(bsx) $(bsx) 1
np = 1
tl = 600
hook = hook/build/libhook.so

include $(shell ap.makesim)

conf: $(hook) b.dat

$(hook): hook
	make -C hook

b.dat: gen
	./$<
