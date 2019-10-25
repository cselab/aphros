all: run

clean:
	rm -f in out main

pipe:
	rm -f in out
	mkfifo in out

run: main
	./main

main: main.cpp
	g++ -std=c++11 -pthread main.cpp -o main

.PHONY: pipe run all clean
