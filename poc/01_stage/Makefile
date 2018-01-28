.PHONY: all build run clean

all: build run

build:
		g++ -fopenmp -g -s -O3 main.cpp -o main
		#clang++ -fopenmp -g -s -O3 main.cpp -o main

run:
	./main

clean:
	rm -v main
