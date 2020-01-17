# Aphros

Solver for incompressible multiphase flows with surface tension.

## Clone

    git clone git@github.com:cselab/aphros.git
    cd aphros
    git config pull.rebase merges

## Documentation

<https://cselab.github.io/hydro/doc/>

generated from `doc/sphinx`

## Build and install

*   Follow instructions from `<ROOT>/deploy/README.md` to
prepare environment and install dependencies.
*   Configure, build and install

     ```
     cd <ROOT>/src

     ./conf                # release
     # ./confdeb           # debug

     make
     # make target=mfer    # to only build ch.mfer
     ```

*   Run tests

     ```
     make test
     ```

*   Install documentation

     ```
     cd <ROOT>/doc/sphinx
     make               # man pages in $CHPREFIX/man
     make web           # html in <ROOT>/doc/sphinx/build/singlehtml
     ```

## Submitting jobs

    ```
    echo NP > np
    echo TL > tl
    ch.submit ARGS  # submit job with NP tasks and time limit TL minutes
    OMP_NUM_THREADS=N ch.submit ARGS  # submit job with NP tasks total and N tasks per node
    ```

Simulation setups in `sim` provide a makefile to submit jobs that run `ch.mfer`

    ```
    echo NP > np
    echo TL > tl
    make submit  # submit job with NP tasks and time limit TL minutes
    ```

## Other repos


### Cubism-hydro

    git clone git@gitlab.ethz.ch:mavt-cse/Cubism-hydro.git
    git clone git@github.com:cselab/Cubism-hydro.git

### hydro

    git clone https://github.com/divfree/hydro.git

### mfer

    git clone git@gitlab.ethz.ch:mavt-cse/mfer.git
