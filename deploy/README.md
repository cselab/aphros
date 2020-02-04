# Setup environment

Requires having `$HOME/bin` in `PATH`.
Prepare the environment and install dependencies:

*   install all packages except for hypre

     ```
     P=prefix make
     ```

    to `prefix` and files `ap-prefix` and `ap-setenv` to `$HOME/bin`
*   set environment

     ```
     . ap-setenv
     ```

*   install hypre

     ```
     cd hypre && make
     ```

*   (optional) install basilisk from `mfer:/ba/deploy`
*   (optional) install basilisk libraries from `sim/ba` by `cmake . && make install`

# Notes

## Prefix

* Read from `$P` if set
* Else from ap-prefix if found in `PATH`
* Else defaults to `prefix`

## Setenv

* `ap-setenv` is written to `$B` (`B=$HOME/bin` in `util`)
  along with `ap-prefix`
* All other units are installed to `$BB` (`BB=$P/bin` in `util`)
* `ap-setenv` adds `$BB` to `PATH` and makes installed units available.
  This prevents the user from running units in an incomplete environment
  (e.g. modules loaded).

## Rules

* Each folder is a package.
* Packages consist of units.
* Unit is an executable or shell script.
* Unit can be executed. Shell scripts can also be sourced.
* _Machine suffix_ is a unique name associated to a particular machine.
* _System unit_ is a unit located in `PATH`.
* Name of a system unit has prefix `ap-` (_standard prefix_).
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
  - need to source `ap-util` to make functions from `util` available.
* Reserved variables:
  - `P`: install path prefix
  - `B`: default `$HOME/bin`
  - `MS`: machine suffix
  - `PRE`: standard prefix

TODO some places assume `PRE=ap-`
