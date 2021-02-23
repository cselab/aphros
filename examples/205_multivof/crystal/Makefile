include par.make

default:
	make cleanrun

include $(shell ap.makesim)

par.conf par.make: gen par.py
	./gen >> out

gen.copy: gen
	cp $< $@

conf: par.conf add.conf gen.copy

par.py:
	touch $@

clean::
	rm -f par.make par.conf bc.dat bg.dat body.dat gen_num.py gen_par.py gen.copy
