# Deploy

* Define the installation prefix. This installs `ap.setenv`.

  ```
  ./install_setenv [--profile PROFILE] PREFIX
  ```

  where `PREFIX` is the installation prefix (e.g. `$HOME/aphros`).

  If specified, `PROFILE` defines a set of wrappers to load modules and submit
  jobs on a cluster. For example, `daint` is for Piz Daint with the
  Slurm workload manager. More profiles can be added to
  directories `scripts/profiles` and `scripts_gen/profiles`.

* Setup the environment. This loads modules and adds the installation prefix to `PATH`.

  ```
  . ap.setenv
  ```

* Install libraries and tools

  ```
  mkdir build
  cd build
  cmake ..
  make -j4
  make install
  ```

  Alternatively, use `ccmake ..` to configure with a dialog.
