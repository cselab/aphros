#make cleanall ; make -j CC=mpic++ bs=32 config=release bgq=1 qpx=1 sequoia=0 hdf=0 hpm=0 precdiv=1 preclevel=1 accurateweno=1 extra=" -D_WENO3_ -D_NOK_"
make cleanall ; make -j CC=mpic++ bs=32 config=release bgq=1 qpx=1 sequoia=0 hdf=0 hpm=0 precdiv=1 preclevel=1 accurateweno=1 extra="-D_BCLABCLOUDABSORB_"
#make cleanall ; make -j CC=mpic++ bs=32 config=debug bgq=1 qpx=1 sequoia=0 hdf=0 hpm=0 precdiv=1 preclevel=1 accurateweno=1 ap=double extra="-D_BCLABCLOUDABSORB_"
