b = 2

ff = $(subst /,,$(dir $(wildcard nx*/out)))

%.dat : %/out
	cd $(basename $@) ; ch.gettraj $(b) > ../$(basename $@).dat

all: $(ff:=.dat)
	gnuplot traj.gp
	gnuplot kedr.gp

.PHONY: all
