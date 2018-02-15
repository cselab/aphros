# Clone

    git clone git@gitlab.ethz.ch:kpetr/Cubism-hydro.git

# Build

## Hypre

    cd lib/build
    P=prefix ./all

## Cubism-hydro

    cd src
    ./conf
    make

# Run
  
    cd sim/sim01
    ./job

# Other repos

## hydro

    git clone https://github.com/divfree/hydro.git 

## mfer

    git clone git@gitlab.ethz.ch:mavt-cse/mfer.git
