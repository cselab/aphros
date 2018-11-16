include $(BASILISK)/Makefile.defs

CFLAGS += -grid=multigrid3D

main: main.c
	qcc -autolink $(CFLAGS) -grid=octree $(LIBS) -lm $< -o $@

run:; mpirun -oversubscribe -np 1 ./main

.PHONY: run
