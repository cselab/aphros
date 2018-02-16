## Rules for tests

* A test consists of 
  - main.cpp
  - CMakeLists.txt
  - data files (if needed)

* All these files are placed in a separate folder

* No single file is aware of all tests

* Build and run all tests with one command

* Build and run an individual test with one command

* Each test should have an function called Simple() 
  containing the minimal client for a feature being tested

* Binary name starts with `t.` to use in gitignore.
  Using a common name (e.g. `main`) not allowed by CMake.
