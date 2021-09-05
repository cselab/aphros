# Test

## Build and run

All tests

    cmake .
    make
    make test

Single test

    cd [directory]
    cmake .
    make 
    make test

## Rules for tests

* A test consists of 
  - main.cpp
  - CMakeLists.txt
  - data files (if needed)

* All these files are placed in a separate directory

* No single file is aware of all tests

* Build and run all tests with one command

* Build and run an individual test with one command

* Each test should define function `Simple()`
  as a minimal example for the feature being tested
  and put it in namespace `simple`.

* Binary name starts with `t.` to use in gitignore.
  Using a common name (e.g. `main`) not allowed by CMake.
