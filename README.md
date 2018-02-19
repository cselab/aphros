# Cubism-hydro

Distributed version of **hydro** based on **Cubism**.

## Clone

    git clone git@gitlab.ethz.ch:kpetr/Cubism-hydro.git

## Build

### Hypre

    cd lib/build
    ./all

### Cubism-hydro

    cd src
    ./conf
    make

    (on daint use `build.daint` instead)

## Run
  
    cd sim/sim01
    ./job

## Other repos

### hydro

    git clone https://github.com/divfree/hydro.git 

### mfer

    git clone git@gitlab.ethz.ch:mavt-cse/mfer.git
