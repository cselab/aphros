# Install lammps serial library

```
git clone git@github.com:lammps/lammps.git
cd lammps/src
make serial -j4 'CCFLAGS = -O3 -g -fPIC'
mkdir -p "$HOME"/.local/lib "$HOME"/.local/include
cp STUBS/libmpi_stubs.a liblammps_serial.a "$HOME"/.local/lib/
cp library.h "$HOME"/.local/include/
```

# Run

LAMMPS configuration is [in.lj]

```
make run
```
