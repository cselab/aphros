include $(BASILISK)/Makefile.defs

main: main.c
	qcc -autolink $(CFLAGS) $(LIBS) -lm $< -o $@

cleandat:
	rm -vf *.vtk
	rm -rvf o
	rm -vf vis/*.{png,mp4}
	rm -vf done main.s.d main.s out traj.dat
	rm -vf Makefile.tests Makefile.deps
