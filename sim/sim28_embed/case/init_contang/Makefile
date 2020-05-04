m = 32 32 32
hook = hook/build/libhook.so

default:
	make cleanall
	make run

include $(shell ap.makesim)

conf: $(hook)

$(hook): hook
	make -C hook


clean::
	rm -vrf sc add.conf
