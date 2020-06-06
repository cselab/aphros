CH ?= "../../ch/doc/sphinx/build/singlehtml/"

copy:
	(cd "$(CH)" && git rev-parse --short HEAD > revision)
	rsync -Pa "$(CH)" .

commit:
	git commit -m "update from devel `cat revision`"

.PHONY: copy commit
