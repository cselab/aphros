omp = 1
bx = 6
px = 1
cubismnc = 1
compress = 0

simmake = sim.makefile

par.make: Makefile pargen
	./pargen $(bx) $(px) $(cubismnc) $(compress)

submit: par.make
	OMP_NUM_THREADS=$(omp) make -f $(simmake) submit

run: par.make
	OMP_NUM_THREADS=$(omp) make -f $(simmake) run

cleanall:
	touch par.make
	make -f $(simmake) cleanall
	rm -f _m _bs _np _tl par.conf par.make

clean:
	rm -f _m _bs _np _tl par.conf par.make
