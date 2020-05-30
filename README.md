# Aphros

Solver for incompressible multiphase flows with surface tension.

## Clone

    git clone git@github.com:cselab/aphros.git
    cd aphros
    git config pull.rebase merges

## Documentation

<https://cselab.github.io/aphros/doc>

generated from `<ROOT>/doc/sphinx`

## Requirements

C++14, cmake, MPI library, hdf5 library with MPI, python, python-numpy.

## Build and install

Follow instructions from `deploy/README.md` to
prepare environment and install dependencies.

Configure, build and install

     cd <ROOT>/src
     make -j4
     # make target=mfer    # to only build ap.mfer

Run tests

     cd <ROOT>/src
     make test

Install documentation

     cd <ROOT>/doc/sphinx
     make               # man pages in $CHPREFIX/man
     make web           # html in <ROOT>/doc/sphinx/build/singlehtml

## Submit job

    echo NP > np
    echo TL > tl
    ap.submit ARGS  # submit job with NP tasks and time limit TL minutes
    OMP_NUM_THREADS=N ap.submit ARGS  # submit job with NP tasks total and N tasks per node

Simulation setups in `sim` provide a makefile to submit jobs that run `ap.mfer`

    echo NP > np
    echo TL > tl
    make submit  # submit job with NP tasks and time limit TL minutes
