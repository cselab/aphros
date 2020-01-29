# Aphros

Solver for incompressible multiphase flows with surface tension.

## Clone

    git clone git@github.com:cselab/aphros.git
    cd aphros
    git config pull.rebase merges

## Migrate to github

    git remote set-url origin git@github.com:cselab/aphros.git

## Documentation

<https://cselab.github.io/hydro/doc/>

generated from `doc/sphinx`

## Build and install

Follow instructions from `<ROOT>/deploy/README.md` to
prepare environment and install dependencies.

Configure, build and install

     cd <ROOT>/src

     ./conf                # release
     # ./confdeb           # debug

     make
     # make target=mfer    # to only build ch.mfer

Run tests

     make test

Install documentation

     cd <ROOT>/doc/sphinx
     make               # man pages in $CHPREFIX/man
     make web           # html in <ROOT>/doc/sphinx/build/singlehtml

## Add copyright notice

Add copyright notice to C/C++ source files found recursively in current
directory (if `copyright` is not found in first 10 lines of the file)

    ap-applycopyright `ap-findsource`

## Submit job

    echo NP > np
    echo TL > tl
    ch.submit ARGS  # submit job with NP tasks and time limit TL minutes
    OMP_NUM_THREADS=N ch.submit ARGS  # submit job with NP tasks total and N tasks per node

Simulation setups in `sim` provide a makefile to submit jobs that run `ch.mfer`

    echo NP > np
    echo TL > tl
    make submit  # submit job with NP tasks and time limit TL minutes

## Other repos


### Cubism-hydro

    git clone git@gitlab.ethz.ch:mavt-cse/Cubism-hydro.git
    git clone git@github.com:cselab/Cubism-hydro.git

### hydro

    git clone https://github.com/divfree/hydro.git

### mfer

    git clone git@gitlab.ethz.ch:mavt-cse/mfer.git
