m = 16 16 16
bs = 8 8 8
hook = hook/build/libhook.so

$(hook): hook
	make -C hook

conf: $(hook)

include $(shell ap.makesim)
