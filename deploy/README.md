# Deploy

* Select profile and create `ap.setenv`

```
./install_setenv PREFIX
(daint) ./install_setenv --profile daint PREFIX
(euler) ./install_setenv --profile euler PREFIX
```

where `PREFIX` is the installation prefix (e.g. `$(HOME)/aphros`)

* Setup environment

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
