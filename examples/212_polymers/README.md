# Coupling with LAMMPS

Install LAMMPS serial library
```
git clone git@github.com:lammps/lammps.git
cd lammps/src
make yes-MOLECULE serial -j4 'CCFLAGS = -O3 -g -fPIC'
mkdir -p "$HOME"/.local/lib "$HOME"/.local/include
cp STUBS/libmpi_stubs.a liblammps_serial.a "$HOME"/.local/lib/
cp library.h "$HOME"/.local/include/
```
LAMMPS configuration is [in.polymers](in.polymers), aphros configuration is in [std.conf](std.conf).

To run the example and generate PPM images,
```
make run
```

Convert to lammps dumps to VTK:
```
awk -f dump2vtk input *.dump
```
