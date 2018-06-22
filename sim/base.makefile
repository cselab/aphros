m ?= 32 32 32
bs ?= 16 16 16
np ?= 1

all: cleanall run

run: conf
	ch.run

submit: conf
	ch.submit

conf: mesh
	ch.base
	ch.aconf

np: 
	echo $(np) > np

mesh: np
	ch.part $(m) $(bs) `cat np` > mesh.conf

clean::
	rm -vf *.{png}
	rm -vf *.{bin,log}
	rm -vf job.id.last job.status arg job.id
	rm -vf mesh.conf base.conf a.conf np

cleandat::
	rm -vf *.pdf
	rm -vf *.{xmf,h5}
	rm -vf *.{vts,pvd}
	rm -vf stat.dat
	rm -vf partit*.csv
	rm -vf s*.vtk
	rm -vf out
	rm -vf lsf.o*
	rm -vf {vx,vy,vz,p,vf}_*.dat

cleanall: clean cleandat

.PHONY: all run submit mesh clean cleandat base conf np
