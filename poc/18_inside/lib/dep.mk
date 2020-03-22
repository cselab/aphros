bbox.o: bbox.c memory.h err.h bbox.h
err.o: err.c err.h
main.o: main.c bbox.h err.h memory.h predicate.h off.h inside.h
memory.o: memory.c
off.o: off.c err.h memory.h off.h inside.h
ply.o: ply.c err.h memory.h ply.h inside.h
predicate.o: predicate.c err.h predicate.h predicate/main.inc
