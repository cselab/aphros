omp = 1
bx = 1
px = 1
cubismnc = 1
compress = 0

simmake = sim.makefile

include $(shell ch.makesim)

par.conf:
	./pargen $(bx) $(px) $(cubismnc) $(compress)

submit: par.conf
	OMP_NUM_THREADS=$(omp) make -C $(simmake) submit

run: conf
	OMP_NUM_THREADS=$(omp) make -C $(simmake) run

clean:
	make -C $(simmake) cleanall
	rm -f _m _bs _np _tl par.conf b.dat par.make
