.POSIX:
#include conf.mk
CXX = mpicxx
M = main
LINK = $(CXX)

all: $M
.o:
	#$(LINK) $< $(LDFLAGS) -L$(CHPREFIX)/lib -L$(HDF5_ROOT)/lib -lhydro -lHYPRE -lhdf5 -o $@
	$(LINK) $< $(LDFLAGS) -L$(CHPREFIX)/lib -L$(HDF5_ROOT)/lib -lhydro -lHYPRE -o $@

.cpp.o:
	$(CXX) $(CFLAGS) -I$(CHPREFIX)/include/hydro $< -c -o $@

.cpp:
%: %.cpp


.PHONY: all clean
clean:; rm -f $M $(M:=.o)
