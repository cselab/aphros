#!/bin/bash
# File       : make.sh
# Date       : Tue Feb 28 12:07:48 2017
# Author     : Fabian Wermelinger
# Description: Code compilation script for node-layer
# Copyright 2017 ETH Zurich. All Rights Reserved.
make -C ../../../CubismApps/MPCFnode/makefiles/ cleanall;
time make -C ../../../CubismApps/MPCFnode/makefiles/ \
    -j CC=mpic++ \
    bs=32 config=release ap=float hdf=1 fftw=0 alphaeps=0.0 \
    qpxemu=1 preclevel=0 microfusion=2 accurateweno=0 \
    extra='-D_NOK_ -D_BCLABCHASIMPLE_ -ldl'
