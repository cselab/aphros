m = 32 32 32
bs = 32 32 32
hook = hook/build/libhook.so

$(hook): hook
	make -C hook

conf: $(hook)
clean::; rm -f *.ppm
include $(shell ap.makesim)
$(hook): hook/hook.cpp
