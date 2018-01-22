CC = gcc

# all: spdp libspdp.a mytest_spdp
all: libspdp.a mytest_spdp

spdp: spdp_10.c
	$(CC) $(extra) -O3 spdp_10.c -o spdp


libspdp.a: spdp_10.c
#	$(CC) $(extra) -O3 -c spdp_10.c
	$(CC) $(extra) -O0 -g -ggdb -c spdp_10.c
	ar rs libspdp.a spdp_10.o
	ranlib libspdp.a

#mytest_spdp: mytest.c
#	$(CC) -fopenmp -O3 -o mytest_spdp mytest.c spdp_10.c -lm

mytest_spdp: mytest.cpp myspdp.h
	$(CC) $(extra) -fopenmp -g -ggdb -O0 -o mytest_spdp mytest.cpp -L. -lspdp -lm
#	$(CC) $(extra) -fopenmp -O3 -o mytest_spdp mytest.cpp -L. -lspdp -lm

clean:
	rm -f spdp mytest_spdp *.o *.a
