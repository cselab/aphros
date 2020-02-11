m = 64 48 32
bs = 8 8 8
np = 12
tl = 10

runomp:
	OMP_NUM_THREADS=12 make run

runmpi:
	OMP_NUM_THREADS=1 make run

include $(shell ap.makesim)

.PHONY: runomp runmpi
