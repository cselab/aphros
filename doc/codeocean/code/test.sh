#!/bin/sh -u

# All tests pass with MPICH but some may fail with OpenMPI.
# This does not affect production runs without Docker.

cd $REPO/src/build/
ctest -j 4 || true

rsync -av $REPO/src/build/Testing /results/
