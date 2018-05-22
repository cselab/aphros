m ?= 64 64 1
bs ?= 32 32 1

all: run

run:
	ch.run

mesh:
	ch.part $(m) $(bs) `cat np` > mesh.conf

clean:
	rm -vf *.{png}
	rm -vf *.pdf
	rm -vf *.{bin,log}
	rm -vf out job.id.last job.status arg

cleandat:
	rm -vf stat.dat
	rm -vf *.dat

cleanall: clean cleandat

.PHONY: run clean cleandat cleanall
