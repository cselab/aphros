m = 32 32 1
np = 1
bs = 32 32 1
hook = hook/build/libhook.so

$(hook): hook
	make -C hook

conf: $(hook)

include $(shell ap.makesim)
