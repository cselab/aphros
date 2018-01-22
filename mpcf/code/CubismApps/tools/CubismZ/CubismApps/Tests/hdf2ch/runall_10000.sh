#!/bin/bash
set -x #echo on

##################
h5file=/home/chatzidp/gitlab/fabdata/data_010000-p.h5
./runone.sh $h5file

##################
h5file=/home/chatzidp/gitlab/fabdata/data_010000-rho.h5
./runone.sh $h5file

##################
h5file=/home/chatzidp/gitlab/fabdata/data_010000-a2.h5
./runone.sh $h5file

##################
h5file=/home/chatzidp/gitlab/fabdata/data_010000-divU.h5
./runone.sh $h5file

##################
h5file=/home/chatzidp/gitlab/fabdata/data_010000-E.h5
./runone.sh $h5file

##################
h5file=/home/chatzidp/gitlab/fabdata/data_010000-Omegax.h5
./runone.sh $h5file

##################
h5file=/home/chatzidp/gitlab/fabdata/data_010000-Ux.h5
./runone.sh $h5file

