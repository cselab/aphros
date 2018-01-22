#!/bin/bash
set -x #echo on

##################
h5file=data-301-p.h5
./runone.sh $h5file

##################
h5file=data-301-g.h5
./runone.sh $h5file
