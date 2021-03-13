include par.make

include $(shell ap.makesim)

par.conf par.make: gen par.py
	./gen >> out

conf: par.conf add.conf

par.py:
	touch $@

clean::
	rm -f par.make par.conf b.dat bc.dat bg.dat body.dat gen_par.py gen_exp.py gen_var.py
