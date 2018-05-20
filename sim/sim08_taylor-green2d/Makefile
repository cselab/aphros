m ?= 64 64 1
bs ?= 32 32 1

all: run

run:
	ch.run

mesh:
	ch.part $(m) $(bs) `cat np` > mesh.conf

clean:
	rm -vf *.{png}
	rm -vf *.{bin,log}
	rm -vf lsf.o*
	rm -vf out job.id.last job.status pid arg

cleandat:
	rm -vf *.pdf
	rm -vf *.{xmf,h5}
	rm -vf *.{vts,pvd}
	rm -vf stat.dat
	rm -vfr sc

cleanall: clean cleandat

.PHONY: run clean cleandat cleanall
