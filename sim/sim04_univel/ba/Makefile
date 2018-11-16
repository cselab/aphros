include $(BASILISK)/Makefile.defs

main: main.c
	#qcc $(CFLAGS) $(BA_LDFLAGS) $< -o $@
	qcc -autolink $(CFLAGS) -lm $(LIBS) $< -o $@

run:; mpirun -oversubscribe -np 1 ./main

.PHONY: run
