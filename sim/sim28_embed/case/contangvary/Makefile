np = 8

default:
	touch add.conf
	make cleanall
	nice make run

2d:
	echo "" > add.conf
	nice make cleanrun m='64 32 1' bs='16 16 1'

3d:
	echo "include 3d.conf" > add.conf
	nice make cleanrun m='64 32 32' bs='16 16 16'

embed_2d:
	echo "include embed.conf" > add.conf
	nice make cleanrun m='64 32 1' bs='16 16 1'

embed_3d:
	echo "include embed_3d.conf" > add.conf
	nice make cleanrun m='64 32 32' bs='16 16 16'

.PHONY: 2d 3d 2d_embed 3d_embed

include $(shell ap.makesim)

conf: add.conf
