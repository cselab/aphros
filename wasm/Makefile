REPO ?= "../../aphros/"
RSYNC = rsync -a -i --update

copy:
	(cd "$(REPO)" && git rev-parse --short HEAD) > revision
	$(RSYNC) $(REPO)/examples/107_wasm/{electrochem,hydro}* .

commit:
	git commit -m "update from devel `cat revision`"

.PHONY: copy commit
