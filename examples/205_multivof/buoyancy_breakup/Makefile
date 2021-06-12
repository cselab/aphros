nx = 32
bsx = 32
m = $(nx) $(nx) 1
bs = $(bsx) $(bsx) 1
np = 1
tl = 600

include $(shell ap.makesim)

conf: b.dat export

export:
	@nx=$(wordlist 1, 1, $m) ; d=nx`printf %03d $$nx` ; \
	echo "mkdir -p $$d && rsync -a {gen,{s,sm}_*.vtk,std.conf,out,out.conf,stat.dat} $$d/" > $@
	@chmod +x $@
	cat $@

.PHONY: export


b.dat: gen
	./$<
