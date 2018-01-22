make cleanall
make -j CC=mpic++ bs=32 config=release bgq=1 qpx=1 sequoia=0 hdf=0 hpm=0 precdiv=1 preclevel=1 accurateweno=1 extra=" -D_WENO3_ "
