# openmp threads per node
omp = 1
# multiplier for number of blocks
kb = 1
# multiplier for number of cores
kp = 1
# block size
bs = 16
cubismnc = 1
compress = 0

simmake = sim.makefile

conf: Makefile pargen
	./pargen $(kb) $(kp) $(bs) $(cubismnc) $(compress)
	OMP_NUM_THREADS=$(omp) make -f $(simmake) conf

submit: conf
	OMP_NUM_THREADS=$(omp) make -f $(simmake) submit

run: conf
	OMP_NUM_THREADS=$(omp) make -f $(simmake) run

cleanall:
	touch par.make
	make -f $(simmake) cleanall
	rm -f _m _bs _np _tl par.conf par.make

clean:
	rm -f _m _bs _np _tl par.conf par.make

.PHONY: conf
