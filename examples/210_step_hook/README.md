# Step hook

Example of customizing the solver loop with a user-defined function.

Function `StepHook()` in `hook/hook.cpp` is called at every step
and takes a pointer to the solver instance `Hydro<M>*`.

Internally, the hook mechanism relies on the `LD_PRELOAD` variable
to load the library `hook/build/libhook.so`
that overrides the default implementation of `StepHook`.

To run the example,

```
make run
```
