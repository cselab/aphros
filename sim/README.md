# Workflow

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
