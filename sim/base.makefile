# mesh size
m ?= 32 32 32
# block size
bs ?= 16 16 16
# number of tasks
np ?= 1
# time limit in minutes
tl ?= 60
# number of threads for openmp
OMP_NUM_THREADS ?= 1

hook ?=
LD_PRELOAD += $(hook)

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
	LD_PRELOAD="$(LD_PRELOAD)" ap.run ap.mfer

submit: conf
	LD_PRELOAD="$(LD_PRELOAD)" ap.submit ap.mfer

kill:
	ap.kill

add.conf:
	touch $@

base.conf:
	ap.base

a.conf:
	ap.aconf

conf: a.conf base.conf mesh.conf add.conf std.conf tl

np:
	echo $(np) > np

tl:
	echo $(tl) > tl

mesh.conf: np
	ap.part $(m) $(bs) `cat np` $(OMP_NUM_THREADS) > mesh.conf

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
	rm -vf sm_*.vtk
	rm -vf bc.vtk
	rm -vf out
	rm -vf lsf.o*
	rm -vf slurm*.out
	rm -vf {vx,vy,vz,p,vf,cl,cls,div,omm}_*.dat
	rm -vf core.*
	rm -vf bc.vtk bc_groups.dat eb.vtk

cleanall: clean cleandat

.PHONY: error cleanrun run submit clean cleandat base conf np tl kill
