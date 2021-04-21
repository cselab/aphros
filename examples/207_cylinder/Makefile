include $(shell ap.makesim)

par.py:
	touch $@

par.conf: par.py gen
	./gen

conf: par.conf

clean::
	rm -vf par.make
