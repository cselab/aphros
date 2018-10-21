m ?= 32 32 32   # mesh size
bs ?= 16 16 16  # block size
np ?= 1         # number of tasks
tl ?= 60        # time limit in minutes

error:
	@echo Error: no target specified. Available targets:
	@echo - cleanrun: cleanall, run
	@echo - run: start in foreground
	@echo - submit: submit job or start in background
	@echo - clean: remove logs, generated conf
	@echo - cleandat: remove output data
	@echo - cleanall: clean, cleandat
	@exit 2


cleanrun: cleanall run

run: conf
	ch.run ch.mfer

submit: conf
	ch.submit ch.mfer

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
	rm -vf mesh.conf base.conf a.conf par.conf np tl

cleandat::
	rm -vf *.pdf
	rm -vf *_*.{xmf,h5}
	rm -vf *_*.vts
	rm -vf p.pvd
	rm -vf stat.dat
	rm -vf partit_*.csv
	rm -vf traj_*.csv
	rm -vf trajsh_*.csv
	rm -vf s_*.vtk
	rm -vf sp_*.vtk
	rm -vf out
	rm -vf lsf.o*
	rm -vf slurm*.out
	rm -vf {vx,vy,vz,p,vf}_*.dat

cleanall: clean cleandat

.PHONY: error cleanrun run submit mesh clean cleandat base conf np tl kill
