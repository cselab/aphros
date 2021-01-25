m = 16 16 16
hook = hook/build/libhook.so

$(hook): hook
	make -C hook

conf: $(hook)

include $(shell ap.makesim)
