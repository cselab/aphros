PUBLICREPO = ../aphros-public

default: help

help:
	@echo 'Error: no target specified. Available targets:'
	@echo '- commit: commit with message from `revision`'
	@echo '- help: display this message'
	@echo
	@exit 2

commit:
	git commit -em "update from devel `head -n1 revision`"

.PHONY: commit default help

