URL = https://github.com/cselab/aphros.git
REPO = aphros

default: clone revision clean

clone:
	git clone --quiet --single-branch --depth 1 $(URL) $(REPO)

revision:
	(cd $(REPO) && git rev-parse HEAD) > $@
	(cd $(REPO) && git log -1 --pretty=%B) >> $@
	(cd $(REPO) && git remote get-url origin) >> $@

clean:
	rm -rf $(REPO)/.git

.PHONY: default clone revision clean
