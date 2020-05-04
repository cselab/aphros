m = 128 128 1
bs = 16 16 1
np = 8
tl = 1440
hook = hook/build/libhook.so

default: conf
.PHONY: default

$(hook): hook
	make -C hook

par.conf:
	nx=`cut -d' ' -f1 <<< "$m"` ; \
	dt=`bc -l <<< "0.5 / $$nx"` ; \
	echo set double dt0 $$dt > par.conf
.PHONY: par.conf

conf: $(hook) par.conf

include $(shell ap.makesim)

clean::
	rm -f omz_*.dat par.conf
