URL = https://github.com/cselab/aphros.git
REPO = aphros

default: clone revision

revision:
	git -C $(REPO) rev-parse HEAD > $@
	git -C $(REPO) log -1 --pretty=%B >> $@
	git -C $(REPO) remote get-url origin >> $@

clone:
	git clone --quiet --single-branch --depth 1 $(URL) $(REPO)
	rm -rf $(REPO)/.git

.PHONY: default clone revision
