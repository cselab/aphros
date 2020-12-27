REPO ?= "../../ch/"
RSYNC = rsync -a -i

copy:
	(cd "$(REPO)" && git rev-parse --short HEAD > revision)
	$(RSYNC) "$(REPO)/doc/sphinx/build/singlehtml/" .
	$(RSYNC) -R $(REPO)/doc/./dashboard/*.{js,html,css} .

commit:
	git commit -m "update from devel `cat revision`"

.PHONY: copy commit
