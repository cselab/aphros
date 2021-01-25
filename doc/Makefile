REPO ?= "../../aphros/"
RSYNC = rsync -a -i --update

copy:
	(cd "$(REPO)" && git rev-parse --short HEAD) > revision
	$(RSYNC) "$(REPO)/doc/sphinx/build/singlehtml/" .

commit:
	git commit -m "update from devel `cat revision`"

.PHONY: copy commit
