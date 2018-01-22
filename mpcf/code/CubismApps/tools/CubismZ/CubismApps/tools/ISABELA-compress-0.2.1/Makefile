SRC_DIR=src 
BSPLINE_DIR=bspline 
WAVELETS_DIR=wavelets
APPS_DIR=apps/example_file
LIB_DIR=lib

all:
	mkdir -p ${LIB_DIR};
	cd ${WAVELETS_DIR}; make;
	cd ${BSPLINE_DIR}; make;
	cd ${SRC_DIR}; make;
	ar rcs ${LIB_DIR}/libisabela.a ${LIB_DIR}/*.o
	rm -f ${LIB_DIR}/*.o
	cd ${APPS_DIR}; make;

clean:
	rm -f ${LIB_DIR}/*.o;
	rm -f ${LIB_DIR}/*.a;
	cd ${WAVELETS_DIR}; make clean;
	cd ${BSPLINE_DIR}; make clean;
	cd ${SRC_DIR}; make clean;
	cd ${APPS_DIR}; make clean;
