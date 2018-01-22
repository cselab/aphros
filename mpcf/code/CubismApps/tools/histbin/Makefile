CC = mpic++
MYFLAGS += -D_SWAP_BYTES_ -I${HOME}/usr/gsl/include 
MYLIBS = -L${HOME}/usr/gsl/lib  -lgsl -lgslcblas

all: printbin bin2txt chunked2ordered printoutliers cleanbin labload

printbin: printbin.cpp
	$(CC) $(MYFLAGS) -o printbin printbin.cpp gnuplot_i.c $(MYLIBS)

bin2txt: bin2txt.cpp
	$(CC) $(MYFLAGS) -o bin2txt bin2txt.cpp

chunked2ordered: chunked2ordered.cpp
	$(CC) $(MYFLAGS) -o chunked2ordered chunked2ordered.cpp

printoutliers: printoutliers.cpp
	$(CC) $(MYFLAGS) -o printoutliers printoutliers.cpp gnuplot_i.c 

cleanbin: cleanbin.cpp
	$(CC) $(MYFLAGS) -o cleanbin cleanbin.cpp

labload: labload.cpp
	$(CC) $(MYFLAGS) -o labload labload.cpp gnuplot_i.c 

clean:
	rm -f printbin bin2txt chunked2ordered printoutliers cleanbin labload
