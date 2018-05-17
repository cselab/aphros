(if not set, default P=prefix)

# Run all packages
    P=prefix ./all

# Run a single package (e.g. hypre)
    cd hypre
    P=prefix ./_run 

`P=prefix` affects only package `setenv` if `ch.setenv` does not exist

# Rules

* Each folder is a package.
* Packages consist of units.
* Unit is an executable or shell script.
* Unit can be executed. Shell scripts can also be sourced.
* _Machine suffix_ is a unique name associated to a particular machine.
* _System unit_ is a unit located in `PATH`.
* Name of a system unit has prefix `.ch` (_standard prefix_).
* System units can be executed or sourced from any folder.
* Unit can call other units by their name. Lookup order:
  - its folder with machine suffix,
  - its folder,
  - `common` with machine suffix,
  - `common`,
  - `PATH` with standard prefix.
* File `util` defines functions available to all units
  - `run()`: lookup and execute,
  - `src()`: lookup and source.
* User can call or source:
  - unit `_run` from containing folder,
  - system units.
* Unit `_run` should source `util`.
* Restrictions apply to system units:
  - can execute or source only system units,
  - need to source `ch.util` to make functions from `util` available.

* Reserved variables:
  - `P`: install path prefix
  - `B`: default `$HOME/bin`
  - `MS`: machine suffix
  - `PRE`: standard prefix
  
TODO some places assume `PRE=ch` 
