nx = 64
bsx = 64
m = $(nx) $(nx) 1
bs = $(bsx) $(bsx) 1
np = 1
tl = 600
hook = hook/build/libhook.so

include $(shell ap.makesim)

conf: $(hook) b.dat export

$(hook): hook
	make -C hook

export:
	d=nx`printf %03d $(nx)` ; \
	echo "mkdir -p $$d && rsync -a {gen,sm_*.vtk,std.conf,out,out.conf,stat.dat} $$d/" > $@
	@chmod +x $@
	cat $@

.PHONY: export

b.dat: gen
	./$<
