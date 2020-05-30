CH ?= "../../ch/doc/sphinx/build/singlehtml/"

copy:
	(cd "$(CH)" && git rev-parse --short HEAD > rev)
	rsync -Pa "$(CH)" .

commit:
	git commit -m "update from hydro revision `cat rev`"

.PHONY: copy commit
