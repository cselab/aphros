include $(shell ap.makesim)

par.py:
	touch $@

par.conf: par.py gen
	./gen

gen.copy: gen
	cp $< $@

conf: par.conf gen.copy

clean::
	rm -vf par.make partinit.csv nucl.csv gen.copy
