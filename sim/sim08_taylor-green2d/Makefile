m ?= 64 64 1
bs ?= 32 32 1

all: run

run: mesh
	ch.run

submit: mesh
	ch.submit

mesh:
	ch.part $(m) $(bs) `cat np` > mesh.conf

clean:
	rm -vf *.{png}
	rm -vf *.{bin,log}
	rm -vf lsf.o*
	rm -vf job.id.last job.status arg

cleandat:
	rm -vf *.pdf
	rm -vf *.{xmf,h5}
	rm -vf *.{vts,pvd}
	rm -vf stat.dat
	rm -vfr sc

cleanall: clean cleandat

.PHONY: all run submit mesh clean cleandat
