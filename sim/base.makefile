m ?= 32 32 32   # mesh size
bs ?= 16 16 16  # block size
np ?= 1         # number of tasks
tl ?= 60        # time limit in minutes

all: cleanall run

run: conf
	ch.run

submit: conf
	ch.submit

kill:
	ch.kill

conf: mesh
	ch.base
	ch.aconf

np: 
	echo $(np) > np

tl: 
	echo $(tl) > tl

mesh: np tl
	ch.part $(m) $(bs) `cat np` > mesh.conf

clean::
	rm -vf *.{png}
	rm -vf *.{bin,log}
	rm -vf job.id.last job.status arg job.id
	rm -vf mesh.conf base.conf a.conf np tl

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

.PHONY: all run submit mesh clean cleandat base conf np tl kill
