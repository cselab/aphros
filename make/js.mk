M4 = m4
M4FLAGS =
.in.js.js:
	@echo creating $@
	@trap 'rm -f $$t; exit 2' 1 2 3 15 && \
	t=/tmp/t.$$$$ && \
	$(M4) $(M4FLAGS) $< > $$t && \
	mv $$t $@
.SUFFIXES: .js .in.js
