m = 32 32 1
bs = 32 32 1
np = 1
tl = 600

#m = 64 64 1
#np = 4

#m = 128 128 1
#np = 8

include $(shell ap.makesim)

conf: b.dat export

export:
	@nx=$(wordlist 1, 1, $m) ; d=nx`printf %03d $$nx` ; \
	echo "mkdir -p $$d && rsync -a {gen,{s,sm}_*.vtk,std.conf} $$d/" > $@
	@chmod +x $@
	cat $@

.PHONY: export


b.dat: gen
	./$<
