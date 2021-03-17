Interactive
[open-source](https://github.com/cselab/aphros/tree/master/examples/107_wasm)
fluid solver [Aphros Explorer](aphros.html) running in the web browser.

# Language

The input configuration is defined using commands
```
[FROMSLIDER] set TYPE NAME VALUE
```
where `TYPE` is `int`, `double`, `vect` or `string`.
Expressions `$VAR_NAME` and `${VAR_NAME}` are expanded to the value
of the variable `VAR_NAME` that must be already defined.
Lines that depend on expansions of variables defined by sliders
need to start with `FROMSLIDER`.

Examples:
```
set int i 0
set double d 1.5
set vect v 1 2 3
set string s text
set string s_m "
  line1
  line2 $s
"
set string sub $i $v $s
FROMSLIDER set vect gravity $g 0
```

# Parameters

## Type `int`

`enable_embed` (0, 1)
: Enable solid bodies using embedded boundaries (cut cells)

`enable_advection` (0, 1)
: Enable advection solver for two-component flows

`hypre_periodic_x`, `hypre_periodic_y` (0, 1)
: Make the domain periodic in the x- and y-directions

`steps_per_frame`
: Number of time steps per rendering frame

`stokes` (0, 1)
: Disable convective fluxes to approach the Stokes limit of no inertia.
  High viscosity or steady state still needed

## Type `double`

`cfl`, `cfla`, `cflst`
: CFL numbers for fluid, advection, and surface tension

`extent`
: Domain size

`mu1`, `mu2`
: Viscosities of components 1 and 2

`rho1`, `rho2`
: Densities of components 1 and 2

`sigma`
: Surface tension

## Type `vect`

`force`
: Body force

`gravity`
: Gravitational acceleration

## Type `string`

`advection_solver` (`vof`, `vofm`)
: Advection solver to use.
Standard volume-of-fluid (`vof`)
and multilayer volume-of-fluid for coalescence prevention (`vofm`)

`bc_path`
: Boundary conditions (see below)

`eb_list_path`
: Primitives composing the solid body (see below)

`init_vf`, list_path`
: Initial volume fraction and list of primitives (see below)

`sliders`
: Sliders (see below)

`visual`
: Fields to visualize (see below)

`print_vars`
: Names of variables to print in the output window


## Other parameters

[Base configuration file](https://github.com/cselab/aphros/blob/master/deploy/scripts/sim_base.conf)
lists other available parameters and their default values.

# Visualization

Parameter `visual` defines a list of scalar fields to visualize.

```
set string visual "
FIELD {
  set vect values VALUES
  set vect colors COLORS
  set vect opacities OPACITIES
}
...
"
```

Valid values of `FIELD`:

`pressure`
: Pressure

`vorticity`
: Vorticity

`embed fraction`
: Fraction of cells occupied by fluid, zero inside solid body

`volume fraction`
: Fraction of cells occupied by second component (gas or liquid)

`velocity magnitude`
: Velocity magnitude

`velocity x`
: Component `x` of velocity

`velocity y`
: Component `y` of velocity

Fields are rendered sequentially, with the colormap and opacities defined by

* field values in `VALUES`,
* flattened RGB color components between 0 and 1 in `COLORS`,
* opacities between 0 and 1 in `OPACITIES`.

Predefined colors that expand to RGB components:

```
set string Cred 1 0.12 0.35
set string Cgreen 0 0.8 0.42
set string Cblue 0 0.6 0.87
set string Cpurple 0.686 0.345 0.729
set string Cyellow 1 0.776 0.118
set string Corange 0.949 0.522 0.133
set string Cwhite 1 1 1
set string Cblack 0 0 0
set string Cgray 0.627 0.694 0.729
```

Example with the `x`-component of velocity overlayed by the volume fraction:

```
set string visual "
  velocity x {
    set vect values -1 0 1
    set vect colors $Cblue $Cwhite $Cred
  }

  volume fraction {
    set vect values 0 1
    set vect colors $Cgreen $Cgreen
    set vect opacities 0 0.5
  }
"
```


# Geometric primitives

Primitives are used to specify initial and boundary conditions.
Valid values of `PRIMITIVE`:

`sphere CX CZ 0 RX [RY=RX]`
: Ellipsoid centered at `CX,CY` with principal semi-axes `RX,RY`

`box CX CZ 0 RX [RY=RX]`
: Box centered at `CX,CY` with half-size `RX,RY`

Examples:

```
sphere 0 0 0 1
sphere 0 1 0 1 2
box 0 0 0 1
box 0 1 0 1 2
```

# Initial volume fraction

Parameter `init_vf` defines the initial volume fraction field.
Valid values of `init_vf`:

`zero`
: zero field

`list`
: initial field defined in `list_path`

Parameter `list_path` defines a list of primitives
```
set string list_path "inline
  PRIMITIVE
  ...
"
```
The initial volume fraction is 1 inside the primitives and 0 elsewhere.

# Solid bodies

Flows in complex geometries around solid bodies are activated by
```
set int enable_embed 1
```
Parameter `eb_list_path` defines a list of primitives that compose the solid body
```
set string eb_list_path "inline
  PRIMITIVE
  ...
"
```

# Boundary conditions

Parameter `bc_path` defines boundary conditions

```
set string bc_path "inline
BCTYPE {
  PRIMITIVE
  ...
}
...
"
```

The boundary conditions are set to `BRTYPE` inside the primitives.
Valid values of `BRTYPE`:

`wall VX VY 0`
: No-slip wall with given velocity

`slipwall`
: Free-slip wall


Example:

```
set string bc_path "inline
wall 0 0 0 {
  box 0 0 0 100
}
wall 1 0 0 {
  box 0 0 0 100 1e-6
  box 0 1 0 100 1e-6
}
"
```

# Sliders

Parameter `sliders` defines a list of sliders

```
set string sliders "
  SLIDER
  ...
"
```

where `SLIDER` describes one slider
```
  VAR MIN MAX [LABEL=VAR]
```
that sets a variable named `VAR` with the range between `MIN` and `MAX`.
Initially, the slider takes the value from the current configuration if
defined and the default `(MIN+MAX)/2` otherwise.
Variables set by sliders are expanded in the input configuration.
Lines that depend on such expansions need to start with `FROMSLIDER`.


Examples:
```
set string sliders "
  rho1 1 10 density 1
  sigma 0 1 surface tension
  g 0 10 gravity
"
FROMSLIDER set vect gravity 0 -$g
```
