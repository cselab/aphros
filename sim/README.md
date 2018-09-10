## List of setups

| id | name                               | ch| ge| 2d| 3d| parameters      | bc               |
|:--:|:-----------------------------------|:-:|:-:|:-:|:-:|:----------------|:-----------------|
| 01 | curvature                          | x | x | x | x | -               | -                |
| 02 | hydrostatic                        | x | x | x | x | Ga              | wall             |
| 03 | cavity                             | x | x | x |   | Re              | wallvel          |
| 04 | univel                             | x | x | x | x | La, Ca          | per, wall        |
| 05 | Mortazavi (from univel)            | x | x | x | x | Re, Ca, x0, R   | per, wall        |
| 06 | Taylor-Green with bubbles          | x | x | x | x | Re, Ca, `b.dat` | per              |
| 07 | growth                             | x |   | x | x | Re, Ca, rate    | per, wall, inlet |


### Boundary conditions

**wallvel**: no-slip with given velocity

**wall**: no-slip with zero velocity

**per**: periodic


# Old

## Workflow

* folders `sim*` contain simulation setups
* setups are machine-independent
* tools to run are installed globally
* tools include:
  - setenv: set environment for build and run
  - run: run locally
  - submit: submit job (or run in background)
  - kill: remove job (or stop background process)
* tool setenv:
  - set environment variables and load modules
  - should be shared with setenv from build
  - supposed to be sourced to current shell (needed for modules)
* tool run:
  - source setenv (though Euler requires loading modules at submit)
  - run executable (normaly, run the job file)
  - wait for executable to finish
* tool submit:
  - source setenv
  - submit job (or run process in background)
  - write job id (or process id) to job.id
* tool kill:
  - remove job (or kill background process)
  - rename to `job.id.last`
