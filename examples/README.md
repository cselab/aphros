# Examples

## Configuration

A simulation directory typically contains the following files:

* `Makefile` with `include $(shell ap.makesim)` which includes the content of [deploy/scripts/sim_base.makefile](../deploy/scripts/sim_base.makefile)
* `std.conf` with commands setting simulation-specific parameters
* `add.conf` optional file with extra parameters

The simulation is launched by `make run` which

* creates the root configuration file `a.conf` which includes files `base.conf`, `mesh.conf`, `std.conf`, and `add.conf`
* creates the base configuration `base.conf` copied from [deploy/scripts/sim_base.conf](../deploy/scripts/sim_base.conf)
* creates `mesh.conf` describing the mesh and its partitioning
* verifies that `std.conf` exists
* creates an empty `add.conf` if it does not exist

The final configuration is the result of executing the files in the following order

* `base.conf` copy of [deploy/scripts/sim_base.conf](../deploy/scripts/sim_base.conf), created automatically
* `mesh.conf` mesh parameters, created automatically
* `std.conf` simulation-specific configuration, provided in the simulation directory
* `add.conf` user-defined overrides, optional file

The default values are listed in [deploy/scripts/sim_base.conf](../deploy/scripts/sim_base.conf).

## Naming conventions

Examples in this directory are numbered as follows:

**000..099**
: Clients for tools without dependencies on numerical solvers.

**100..199**
: Standalone programs that use numerical solvers from Aphros.

**200..299**
: Simulation setups for the full solver (executable `ap.mfer`).

