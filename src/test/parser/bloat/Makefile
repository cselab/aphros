APHROS = -I${APHROS_PREFIX}/include/ -L${APHROS_PREFIX}/lib -laphros 

main: main.cpp
	time g++ -O3 $(APHROS) $< _bloat*.cpp -o $@

clean:
	rm -f _bloat*.cpp main.cpp all.h

gen:
	./gen.sh

.PHONY: clean gen
