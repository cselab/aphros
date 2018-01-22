#AM_CFLAGS = -I./include -I../zlib
#LDFLAGS=-fPIC -shared

AUTOMAKE_OPTIONS=foreign
include_HEADERS=include/CompressElement.h include/DynamicDoubleArray.h include/rw.h\
		include/dictionary.h include/DynamicFloatArray.h include/VarSet.h include/sz.h include/Huffman.h\
		include/DynamicByteArray.h include/DynamicIntArray.h include/iniparser.h include/TightDataPointStorageD.h include/TightDataPointStorageF.h sz.mod rw.mod

lib_LTLIBRARIES=libsz.la
libsz_la_CFLAGS=-g -I./include -I../zlib/
libsz_la_LDFLAGS = -version-info  1:4:0
libsz_la_LIDADD=../zlib/.libs/libzlib.a
libsz_la_SOURCES=src/ByteToolkit.c src/dataCompression.c src/DynamicIntArray.c src/iniparser.c src/szf.c \
		src/CompressElement.c src/DynamicByteArray.c src/rw.c src/TightDataPointStorageD.c src/TightDataPointStorageF.c \
		src/conf.c src/DynamicDoubleArray.c src/rwf.c src/TypeManager.c \
		src/dictionary.c src/DynamicFloatArray.c src/VarSet.c src/test_zlib.c src/Huffman.c src/sz_float.c src/sz_double.c src/sz.c \
		src/sz_interface.F90 src/rw_interface.F90
