np = 1

MAKE=make -f sim.makefile

export m bs

runpost:
	$(MAKE) cleanrun
	$(MAKE) post

post:
	ap.getcol c2.x > c2x

2d:
	echo "include 2d.conf" > add.conf
	$(MAKE) runpost m='16 16 1' bs='16 16 1'

3d:
	echo "include 3d.conf" > add.conf
	$(MAKE) runpost m='16 16 16' bs='16 16 16'

embed_2d:
	echo "include embed_2d.conf" > add.conf
	$(MAKE) runpost m='16 16 1' bs='16 16 1'

embed_3d:
	echo "include embed_3d.conf" > add.conf
	$(MAKE) runpost m='16 16 16' bs='16 16 16'

.PHONY: 2d 3d embed_2d embed_3d runpost post

include $(shell ap.makesim)
