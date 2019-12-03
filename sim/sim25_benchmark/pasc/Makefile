omp = 1
bx = 1
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

clean:
	touch par.make
	make -f $(simmake) cleanall
	rm -f _m _bs _np _tl par.conf par.make
