#!/bin/bash
set -x #echo on

##################
h5file=/home/chatzidp/gitlab/fabdata/data_005000-p.h5
./runone.sh $h5file

##################
h5file=/home/chatzidp/gitlab/fabdata/data_005000-rho.h5
./runone.sh $h5file

##################
h5file=/home/chatzidp/gitlab/fabdata/data_005000-a2.h5
./runone.sh $h5file

##################
h5file=/home/chatzidp/gitlab/fabdata/data_005000-divU.h5
./runone.sh $h5file

##################
h5file=/home/chatzidp/gitlab/fabdata/data_005000-E.h5
./runone.sh $h5file

##################
h5file=/home/chatzidp/gitlab/fabdata/data_005000-Omegax.h5
./runone.sh $h5file

##################
h5file=/home/chatzidp/gitlab/fabdata/data_005000-Ux.h5
./runone.sh $h5file
