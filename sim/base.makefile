m ?= 32 32 32
bs ?= 16 16 16

all: cleanall run

run: conf
	ch.run

submit: conf
	ch.submit

conf: mesh
	ch.base
	ch.aconf

mesh:
	ch.part $(m) $(bs) `cat np` > mesh.conf

clean::
	rm -vf *.{png}
	rm -vf *.{bin,log}
	rm -vf job.id.last job.status arg
	rm -vf mesh.conf base.conf a.conf

cleandat::
	rm -vf *.pdf
	rm -vf *.{xmf,h5}
	rm -vf *.{vts,pvd}
	rm -vf stat.dat
	rm -vf partit*.csv
	rm -vf s*.vtk
	rm -vf out
	rm -vf lsf.o*

cleanall: clean cleandat

.PHONY: all run submit mesh clean cleandat base conf
