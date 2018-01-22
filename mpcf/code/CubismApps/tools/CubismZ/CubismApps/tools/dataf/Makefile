CC = g++
CPPFLAGS = -O3
ap ?= float

ifeq "$(ap)" "float"
	CPPFLAGS += -D_FLOAT_PRECISION_
endif

all: printdat writedat compdat compdat2 printdat_stride


printdat: printdat.cpp
printdat_stride: printdat_stride.cpp
writedat: writedat.cpp
compdat: compdat.cpp
compdat2: compdat2.cpp

clean:
	rm -f printdat writedat compdat compdat2 printdat_stride
