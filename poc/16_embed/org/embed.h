/**
# Embedded boundaries

Boundaries of general shape can be described using an integral
(i.e. finite volume) formulation which takes into account the volume
and area fractions of intersection of the embedded boundary with the
Cartesian mesh. 

We will need to deal with volume fractions. Interpolations (for
Dirichlet boundary conditions) assume a 5x5 stencil. */

#include "fractions.h"
#define BGHOSTS 2
#define EMBED 1

/**
The volume and area fractions are stored in these fields. */

scalar cs[];
face vector fs[];

/**
Embedded boundary operators specific to trees are defined in this
file. */

#if TREE
# include "embed-tree.h"
#endif

/**
## Operator overloading

Several standard operators, defined in [common.h]() need to be tuned
to take into account the embedded fractions. 

The *SEPS* constant is used to avoid division by zero. */

#undef SEPS
#define SEPS 1e-30

/**
When combining third-order Dirichlet conditions and [approximate
projections](navier-stokes/centered.h#approximate-projection),
instabilities may occur due to feedback between the pressure mode and
the velocity, amplified by the third-order derivative. This can be
stabilised using the weighted average below when computing face
velocities. The corresponding test case is [test/uf.c](). Note that if
only second-order Dirichlet fluxes are used, simple averaging is
stable. */

#define cs_avg(a,i,j,k)							\
  ((a[i,j,k]*(1.5 + cs[i,j,k]) + a[i-1,j,k]*(1.5 + cs[i-1,j,k]))/	\
   (cs[i,j,k] + cs[i-1,j,k] + 3.))

/**
Face gradients and face values, computed from cell-centered values
must be tuned to take into account the area fractions of the embedded
boundary. This follows the procedure described in [Johansen and
Colella, 1998](#johansen1998), figure 3 and equation (16) in
particular. Note that this is only used when using second-order
fluxes. */

#if dimension == 2
#define face_condition(fs, cs)						\
  (fs.x[i,j] > 0.5 && fs.y[i,j + (j < 0)] && fs.y[i-1,j + (j < 0)] &&	\
   cs[i,j] && cs[i-1,j])

foreach_dimension()
static inline double embed_face_gradient_x (Point point, scalar a, int i)
{
  int j = sign(fs.x[i,1] - fs.x[i,-1]);
  assert (cs[i] && cs[i-1]);
  if (face_condition (fs, cs))
    return ((1. + fs.x[i])*(a[i] - a[i-1]) +
	    (1. - fs.x[i])*(a[i,j] - a[i-1,j]))/(2.*Delta);
  return (a[i] - a[i-1])/Delta;
}

foreach_dimension()
static inline double embed_face_value_x (Point point, scalar a, int i)
{
  int j = sign(fs.x[i,1] - fs.x[i,-1]);
  return face_condition (fs, cs) ?
    ((1. + fs.x[i])*cs_avg(a,i,0,0) + (1. - fs.x[i])*cs_avg(a,i,j,0))/2. :
    cs_avg(a,i,0,0);
}

/**
The generalisation to 3D is a bit more complicated. See Fig. 1 of
[Schwartz et al, 2006](#schwartz2006). */

#else // dimension == 3
foreach_dimension()
static inline coord embed_face_barycentre_z (Point point, int i)
{
  // Young's normal calculation
  coord n1 = {0};
  double nn = 0.;
  scalar f = fs.z;
  foreach_dimension(2) {
    n1.x = (f[-1,-1,i] + 2.*f[-1,0,i] + f[-1,1,i] -
	    f[+1,-1,i] - 2.*f[+1,0,i] - f[+1,1,i]);
    nn += fabs(n1.x);
  }
  if (!nn)
    return (coord){0.,0.,0.};
  foreach_dimension(2)
    n1.x /= nn;
  // Position `p` of the face barycentre
  coord n, p1, p;
  ((double *)&n)[0] = n1.x, ((double *)&n)[1] = n1.y;
  double alpha = line_alpha (f[0,0,i], n);
  line_center (n, alpha, f[0,0,i], &p1);
  p.x = ((double *)&p1)[0], p.y = ((double *)&p1)[1], p.z = 0.;
  return p;
}

/**
The macro below defines the condition which must be verified when
applying the second-order flux interpolation and/or face
averaging. The goal is to avoid using values from areas of the mesh
which are not topologically connected. Doing so would couple
disconnected subproblems (for example when solving a Laplacian) which
would most probably lead to lack of convergence. This is very
important for robustness when dealing with complex boundaries. */

#define face_condition(fs, cs)						\
  (fs.x[i,j,k] > 0.5 && (fs.x[i,j,0] > 0.5 || fs.x[i,0,k] > 0.5) &&	\
   fs.y[i,j + (j < 0),0] && fs.y[i-1,j + (j < 0),0] &&			\
   fs.y[i,j + (j < 0),k] && fs.y[i-1,j + (j < 0),k] &&			\
   fs.z[i,0,k + (k < 0)] && fs.z[i-1,0,k + (k < 0)] &&			\
   fs.z[i,j,k + (k < 0)] && fs.z[i-1,j,k + (k < 0)] &&			\
   cs[i-1,j,0] && cs[i-1,0,k] && cs[i-1,j,k] &&				\
   cs[i,j,0] && cs[i,0,k] && cs[i,j,k])

foreach_dimension()
static inline double embed_face_gradient_x (Point point, scalar a, int i)
{
  assert (cs[i] && cs[i-1]);
  coord p = embed_face_barycentre_x (point, i);
  // Bilinear interpolation of the gradient (see Fig. 1 of Schwartz et al., 2006)
  int j = sign(p.y), k = sign(p.z);
  if (face_condition(fs, cs)) {
    p.y = fabs(p.y), p.z = fabs(p.z);
    return (((a[i,0,0] - a[i-1,0,0])*(1. - p.y) +
	     (a[i,j,0] - a[i-1,j,0])*p.y)*(1. - p.z) + 
	    ((a[i,0,k] - a[i-1,0,k])*(1. - p.y) +
	     (a[i,j,k] - a[i-1,j,k])*p.y)*p.z)/Delta;
  }
  return (a[i] - a[i-1])/Delta;
}

foreach_dimension()
static inline double embed_face_value_x (Point point, scalar a, int i)
{
  coord p = embed_face_barycentre_x (point, i);
  // Bilinear interpolation
  int j = sign(p.y), k = sign(p.z);
  if (face_condition(fs, cs)) {
    p.y = fabs(p.y), p.z = fabs(p.z);
    return ((cs_avg(a,i,0,0)*(1. - p.y) + cs_avg(a,i,j,0)*p.y)*(1. - p.z) + 
	    (cs_avg(a,i,0,k)*(1. - p.y) + cs_avg(a,i,j,k)*p.y)*p.z);
  }
  return cs_avg(a,i,0,0);
}
#endif // dimension == 3

/**
We use the functions above to redefine the face gradient macros. Note
that the second-order face gradients and averaging are used only if
the corresponding scalar attribute below (`third` because of
third-order accuracy when using Dirichlet conditions, see
[test/neumann.c]) is set to `true`. The default is `false`. */

attribute {
  bool third;
}

#undef face_gradient_x
#define face_gradient_x(a,i)					\
  (a.third && fs.x[i] < 1. && fs.x[i] > 0. ?			\
   embed_face_gradient_x (point, a, i) :			\
   (a[i] - a[i-1])/Delta)

#undef face_gradient_y
#define face_gradient_y(a,i)					\
  (a.third && fs.y[0,i] < 1. && fs.y[0,i] > 0. ?		\
   embed_face_gradient_y (point, a, i) :			\
   (a[0,i] - a[0,i-1])/Delta)

#undef face_gradient_z
#define face_gradient_z(a,i)					\
  (a.third && fs.z[0,0,i] < 1. && fs.z[0,0,i] > 0. ?		\
   embed_face_gradient_z (point, a, i) :			\
   (a[0,0,i] - a[0,0,i-1])/Delta)

#undef face_value
#define face_value(a,i)							\
  (a.third && fs.x[i] < 1. && fs.x[i] > 0. ?				\
   embed_face_value_x (point, a, i) :					\
   cs_avg(a,i,0,0))

/**
The centered gradient must not use values of fields entirely contained
within the embedded boundary (for which *cs* is zero). */

#undef center_gradient
#define center_gradient(a) (fs.x[] && fs.x[1] ? (a[1] - a[-1])/(2.*Delta) : \
			    fs.x[1] ? (a[1] - a[])/Delta :		    \
			    fs.x[]  ? (a[] - a[-1])/Delta : 0.)

/**
## Utility functions for the geometry of embedded boundaries

For a cell containing a fragment of embedded boundary (i.e. for which
$0 < cs < 1$), *embed_geometry()* returns the area of the fragment,
the relative position *p* of the barycenter of the fragment and the
boundary normal *n*. */

static inline
double embed_geometry (Point point, coord * p, coord * n)
{
  *n = facet_normal (point, cs, fs);
  double alpha = plane_alpha (cs[], *n);
  double area = plane_area_center (*n, alpha, p);
  normalize (n);
  return area;
}

/**
This function and the macro below shift the position $(x1,y1,z1)$ to
the position of the barycenter of the embedded fragment.  */

static inline
double embed_area_center (Point point, double * x1, double * y1, double * z1)
{
  double area = 0.;
  if (cs[] > 0. && cs[] < 1.) {
    coord n, p;
    area = embed_geometry (point, &p, &n);
    *x1 += p.x*Delta, *y1 += p.y*Delta, *z1 += p.z*Delta;
  }
  return area;
}

#define embed_pos() embed_area_center (point, &x, &y, &z)

/**
This function returns the value of field *s* interpolated linearly at
the barycenter *p* of the fragment of embedded boundary contained
within the cell. */

double embed_interpolate (Point point, scalar s, coord p)
{
  assert (dimension == 2);
  int i = sign(p.x), j = sign(p.y);
  if (cs[i] && cs[0,j] && cs[i,j])
    // bilinear interpolation when all neighbors are defined
    return ((s[]*(1. - fabs(p.x)) + s[i]*fabs(p.x))*(1. - fabs(p.y)) + 
	    (s[0,j]*(1. - fabs(p.x)) + s[i,j]*fabs(p.x))*fabs(p.y));
  else {
    // linear interpolation with gradients biased toward the
    // cells which are defined
    double val = s[];
    foreach_dimension() {
      int i = sign(p.x);
      if (cs[i])
	val += fabs(p.x)*(s[i] - s[]);
      else if (cs[-i])
	val += fabs(p.x)*(s[] - s[-i]);
    }
    return val;
  }
}

/**
This function "removes" (by setting their volume fraction to zero)
cells which have inconsistent volume/surface fractions. This is
important to guarantee the robustness of the solution for complex (and
under-resolved) boundaries. The functions returns the number of cells
removed. */

struct Cleanup {
  scalar c;
  face vector s;
  double smin;   // minimum surface fraction (optional)
  bool opposite; // whether to eliminate 'thin tubes' (optional)
};

trace
int fractions_cleanup (struct Cleanup p)
{
  scalar c = p.c;
  face vector s = p.s;
  
  /**
  Since both surface and volume fractions are altered, iterations are
  needed. This reflects the fact that changes are coupled spatially
  through the topology of the domain: for examples, long, unresolved
  "filaments" may need many iterations to be fully removed. */
  
  int changed = 1, schanged = 0;
  for (int i = 0; i < 100 && changed; i++) {

    /**
    Face fractions of empty cells must be zero. */
   
    foreach_face()
      if (s.x[] && ((!c[] || !c[-1]) || s.x[] < p.smin))
	s.x[] = 0.;
    boundary ((scalar *){s});

    changed = 0;
    foreach(reduction(+:changed))
      if (c[] > 0. && c[] < 1.) {
	int n = 0;
	foreach_dimension() {
	  for (int i = 0; i <= 1; i++)
	    if (s.x[i] > 0.)
	      n++;

	  /**
	  If opposite surface fractions are zero (and volume fraction
	  is non-zero), then we are dealing with a thin "tube", which
	  we just remove because it can sometimes lead to
	  non-convergence when
	  [projecting](navier-stokes/centered.h#approximate-projection)
	  the velocity field. */

	  if (p.opposite && s.x[] == 0. && s.x[1] == 0.)
	    c[] = 0., changed++;
	}

	/**
	The number of "non-empty" faces (i.e. faces which have a
	surface fraction larger than epsilon) cannot be smaller than
	the dimension (the limiting cases correspond to a triangle in
	2D and a tetrahedron in 3D). */
	
	if (n < dimension)
	  c[] = 0., changed++;
      }
    boundary ({c});

    schanged += changed;
  }
  assert (!changed); // not converged yet...
  return schanged;
}
  
/**
## Dirichlet boundary condition

This function returns the gradient of scalar *s*, normal to the
embedded boundary defined by *cs*, of unit normal vector *n*
(normalised using the Euclidean norm, not the box norm) and of
centroid *p*. The Dirichlet boundary condition *bc* is imposed on the
embedded boundary.

The calculation follows [Johansen and Colella, 1998](#johansen1998)
and is summarised in the figure below (see also Figure 4 of Johansen
and Colella and Figure 2 of [Schwartz et al, 2006](#schwartz2006) for
the 3D implementation).

![Third-order normal gradient scheme](figures/dirichlet_gradient.svg) 

For degenerate cases, a non-zero value of *coef* is returned and
`coef*s[]` must be added to the value returned to obtain the gradient. */

#define quadratic(x,a1,a2,a3) \
  (((a1)*((x) - 1.) + (a3)*((x) + 1.))*(x)/2. - (a2)*((x) - 1.)*((x) + 1.))

foreach_dimension()
static inline double dirichlet_gradient_x (Point point, scalar s, scalar cs,
					   coord n, coord p, double bc,
					   double * coef)
{
  foreach_dimension()
    n.x = - n.x;
  double d[2], v[2] = {nodata,nodata};
  bool defined = true;
  foreach_dimension()
    if (defined && !fs.x[(n.x > 0.)])
      defined = false;
  if (defined)
    for (int l = 0; l <= 1; l++) {
      int i = (l + 1)*sign(n.x);
      d[l] = (i - p.x)/n.x;
      double y1 = p.y + d[l]*n.y;
      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
      y1 -= j;
#if dimension == 2
      if (fs.x[i + (i < 0),j] && fs.y[i,j] && fs.y[i,j+1] &&
	  cs[i,j-1] && cs[i,j] && cs[i,j+1])
	v[l] = quadratic (y1, (s[i,j-1]), (s[i,j]), (s[i,j+1]));
#else // dimension == 3
      double z = p.z + d[l]*n.z;
      int k = z > 0.5 ? 1 : z < -0.5 ? -1 : 0;
      z -= k;
      bool defined = fs.x[i + (i < 0),j,k];
      for (int m = -1; m <= 1 && defined; m++)
	if (!fs.y[i,j,k+m] || !fs.y[i,j+1,k+m] ||
	    !fs.z[i,j+m,k] || !fs.z[i,j+m,k+1] ||
	    !cs[i,j+m,k-1] || !cs[i,j+m,k] || !cs[i,j+m,k+1])
	  defined = false;
      if (defined)
	// bi-quadratic interpolation
	v[l] =
	  quadratic (z,
		     quadratic (y1,
				(s[i,j-1,k-1]), (s[i,j,k-1]), (s[i,j+1,k-1])),
		     quadratic (y1,
				(s[i,j-1,k]),   (s[i,j,k]),   (s[i,j+1,k])),
		     quadratic (y1,
				(s[i,j-1,k+1]), (s[i,j,k+1]), (s[i,j+1,k+1])));
#endif // dimension == 3
      else
	break;
    }
  if (v[0] == nodata) {

    /**
    This is a degenerate case, we use the boundary value and the
    cell-center value to define the gradient. */
	
    d[0] = max(1e-3, fabs(p.x/n.x));
    *coef = - 1./(d[0]*Delta);
    return bc/(d[0]*Delta);
  }

  /**
  For non-degenerate cases, the gradient is obtained using either
  second- or third-order estimates. */
  
  *coef = 0.;
  if (v[1] != nodata) // third-order gradient
    return (d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta);
  return (bc - v[0])/(d[0]*Delta); // second-order gradient
}

double dirichlet_gradient (Point point, scalar s, scalar cs,
			   coord n, coord p, double bc, double * coef)
{
#if dimension == 2
  foreach_dimension()
    if (fabs(n.x) >= fabs(n.y))
      return dirichlet_gradient_x (point, s, cs, n, p, bc, coef);
#else // dimension == 3
  if (fabs(n.x) >= fabs(n.y)) {
    if (fabs(n.x) >= fabs(n.z))
      return dirichlet_gradient_x (point, s, cs, n, p, bc, coef);
  }
  else if (fabs(n.y) >= fabs(n.z))
    return dirichlet_gradient_y (point, s, cs, n, p, bc, coef);
  return dirichlet_gradient_z (point, s, cs, n, p, bc, coef);
#endif // dimension == 3
  return nodata;
}

bid embed;

/**
## Surface force and vorticity

We first define a function which computes
$\mathbf{\nabla}\mathbf{u}\cdot\mathbf{n}$ while taking the boundary
conditions on the embedded surface into account. */

static inline
coord embed_gradient (Point point, vector u, coord p, coord n)
{
  coord dudn;
  foreach_dimension() {
    bool dirichlet;
    double vb = u.x.boundary[embed] (point, point, u.x, &dirichlet);
    if (dirichlet) {
      double val;
      dudn.x = dirichlet_gradient (point, u.x, cs, n, p, vb, &val);
    }
    else // Neumann
      dudn.x = vb;
    if (dudn.x == nodata)
      dudn.x = 0.;
  }
  return dudn;
}

/**
The force exerted by the fluid on the solid can be written
$$
\mathbf{F}_{\Gamma} = - \int_{\partial \Gamma} ( - p\mathbf{I} +
2 \mu \mathbf{D}) \cdot \mathbf{n}d \partial \Gamma
$$
with $\Gamma$ the solid boundary. It can be further decomposed into a
pressure (i.e. "form") drag
$$
\mathbf{F}_p = \int_{\partial \Gamma} p \mathbf{n}d \partial \Gamma
$$
and a viscous drag
$$
\mathbf{F}_{\mu} = - \int_{\partial \Gamma} 
2 \mu \mathbf{D} \cdot \mathbf{n}d \partial \Gamma
$$
These two vectors are computed by the *embed_force()* function.
*/

trace
void embed_force (scalar p, vector u, face vector mu, coord * Fp, coord * Fmu)
{
  // fixme: this could be simplified considerably if reduction worked on vectors
  double Fpx = 0., Fpy = 0., Fmux = 0., Fmuy = 0.;
  foreach (reduction(+:Fpx) reduction(+:Fpy)
	   reduction(+:Fmux) reduction(+:Fmuy))
    if (cs[] > 0. && cs[] < 1.) {

      /**
      To compute the pressure force, we first get the coordinates of
      the barycentre of the embedded fragment, its area and normal,
      and then interpolate the pressure field on the surface. */
      
      coord n, b;
      double area = embed_geometry (point, &b, &n);
      area *= pow (Delta, dimension - 1);
      double Fn = area*embed_interpolate (point, p, b);
      Fpx += Fn*n.x;
      Fpy += Fn*n.y;

      /**
      To compute the viscous force, we first need to retrieve the
      local value of the viscosity (ideally at the barycentre of the
      embedded fragment). This is not completely trivial since it is
      defined on the faces of the cell. We use a
      surface-fraction-weighted average value. */
      
      if (constant(mu.x) != 0.) {
	double mua = 0., fa = 0.;
	foreach_dimension() {
	  mua += mu.x[] + mu.x[1];
	  fa  += fs.x[] + fs.x[1];
	}
	mua /= fa;

	/**
	To compute the viscous force, we need to take into account the
	(Dirichlet) boundary conditions for the velocity on the
	surface. We only know how to do this when computing the normal
	gradient $\mathbf{\nabla}\mathbf{u}\cdot\mathbf{n}$ using the
	[embed_gradient()](#embed_gradient) function. We thus
	need to re-express the viscous force using only normal
	derivatives of the velocity field.
	
	If we assume that $\mathbf{u}$ is constant on the boundary, then
	$$
	\mathbf{{\nabla}} \mathbf{u} \cdot \mathbf{t}= \mathbf{0}
	$$
	with $\mathbf{t}$ the unit tangent vector to the boundary. We
	thus have the relations
	$$
	\mathbf{{\nabla}} \mathbf{u} = \left( \mathbf{{\nabla}} \mathbf{u}
	\cdot \mathbf{n} \right) \mathbf{n} + \left( \mathbf{{\nabla}}
	\mathbf{u} \cdot \mathbf{t} \right) \mathbf{t} = \left(
	\mathbf{{\nabla}} \mathbf{u} \cdot \mathbf{n} \right) \mathbf{n}
	$$
	$$
	\mathbf{D}= \frac{1}{2}  \left( \mathbf{{\nabla}} \mathbf{u} +
	\mathbf{{\nabla}}^T \mathbf{u} \right) = \frac{1}{2} 
	\left(\begin{array}{cc}
	2 \left( \mathbf{{\nabla}} u \cdot \mathbf{n} \right) n_x & \left(
	\mathbf{{\nabla}} u \cdot \mathbf{n} \right) n_y + \left(
	\mathbf{{\nabla}} v \cdot \mathbf{n} \right) n_x\\
	\left( \mathbf{{\nabla}} u \cdot \mathbf{n} \right) n_y + \left(
	\mathbf{{\nabla}} v \cdot \mathbf{n} \right) n_x & 2 \left(
	\mathbf{{\nabla}} v \cdot \mathbf{n} \right) n_y
	\end{array}\right)
	$$
	$$
	\mathbf{F}_{\mu} = - \int_{\Gamma} \left(\begin{array}{c}
	\left[2 \mu \left( \mathbf{{\nabla}} u \cdot \mathbf{n} \right)
	n_x \right] n_x + \mu \left[ \left( \mathbf{{\nabla}} u \cdot \mathbf{n}
	\right) n_y + \left( \mathbf{{\nabla}} v \cdot \mathbf{n} \right) n_x
	\right] n_y\\
	\left[2 \mu \left( \mathbf{{\nabla}} v \cdot \mathbf{n} \right)
	n_y \right] n_y + \mu \left[ \left( \mathbf{{\nabla}} u \cdot \mathbf{n}
	\right) n_y + \left( \mathbf{{\nabla}} v \cdot \mathbf{n} \right) n_x
	\right] n_x
	\end{array}\right)
	$$
	$$
	\mathbf{F}_{\mu} = - \int_{\Gamma} \left(\begin{array}{c}
	\mu \left[ \left( \mathbf{{\nabla}} u \cdot \mathbf{n} \right) 
	(n^2_x + 1) + \left( \mathbf{{\nabla}} v \cdot \mathbf{n} \right) n_x
	n_y \right]\\
	\mu \left[ \left( \mathbf{{\nabla}} v \cdot \mathbf{n} \right) 
	(n^2_y + 1) + \left( \mathbf{{\nabla}} u \cdot \mathbf{n} \right) n_x
	n_y \right]
	\end{array}\right)
	$$
	*/

	assert (dimension == 2);
	coord dudn = embed_gradient (point, u, b, n);
	Fmux -= area*mua*(dudn.x*(sq(n.x) + 1.) + dudn.y*n.x*n.y);
	Fmuy -= area*mua*(dudn.y*(sq(n.y) + 1.) + dudn.x*n.x*n.y);
      }
    }

  Fp->x = Fpx; Fp->y = Fpy; 
  Fmu->x = Fmux; Fmu->y = Fmuy; 
}

/**
In two dimensions, *embed_vorticity()* returns the vorticity of
velocity field *u*, on the surface of the embedded boundary contained
in the cell. *p* is the relative position of the barycentre of the
embedded fragment and *n* its normal. */

#if dimension == 2
double embed_vorticity (Point point, vector u, coord p, coord n)
{
  /**
  We compute $\mathbf{{\nabla}}\mathbf{u}\cdot\mathbf{n}$, taking
  the boundary conditions into account. */
    
  coord dudn = embed_gradient (point, u, p, n);

  /**
  The vorticity is then obtained using the relations
  $$
  \omega = \partial_x v - \partial_y u = 
  \left( \mathbf{{\nabla}} v \cdot \mathbf{n} \right) n_x - 
  \left( \mathbf{{\nabla}} u \cdot \mathbf{n} \right) n_y
  $$
  */
    
  return dudn.y*n.x - dudn.x*n.y;
}
#endif // dimension == 2

/**
## Flux through the embedded boundary

This function computes the flux through the embedded boundary contained 
within a cell
$$
\int_b \mu \nabla s\cdot\mathbf{n} db
$$
with $db$ the elementary boundary surface and $\mathbf{n}$ the embedded
boundary (outward-pointing) normal.

Boundary conditions for *s* are taken into account.

The result is returned in *val*. 

For degenerate cases, the value returned by the function must be
multiplied by `s[]` and added to *val*. */

double embed_flux (Point point, scalar s, face vector mu, double * val)
{

  /**
  If the cell does not contain a fragment of embedded boundary, the
  flux is zero. */
  
  *val = 0.;
  if (cs[] >= 1. || cs[] <= 0.)
    return 0.;

  /**
  If the boundary condition is homogeneous Neumann, the flux is
  zero. */

  bool dirichlet;
  double grad = s.boundary[embed] (point, point, s, &dirichlet);
  if (!grad && !dirichlet)
    return 0.;

  /**
  We compute the normal, area and barycenter of the fragment of embedded
  boundary contained within the cell. */

  coord n = facet_normal (point, cs, fs), p;
  double alpha = plane_alpha (cs[], n);
  double area = plane_area_center (n, alpha, &p);

  /**
  If the boundary condition is Dirichlet, we need to compute the
  normal gradient. */

  double coef = 0.;
  if (dirichlet) {
    normalize (&n);
    grad = dirichlet_gradient (point, s, cs, n, p, grad, &coef);
  }

  /**
  We retrieve the (average) value of $\mu$ without the metric. */
  
  double mua = 0., fa = 0.;
  foreach_dimension() {
    mua += mu.x[] + mu.x[1];
    fa  += fs.x[] + fs.x[1];
  }
  *val = - mua/(fa + SEPS)*grad*area/Delta;
  return - mua/(fa + SEPS)*coef*area/Delta;;
}

/**
For ease of use, we redefine the Neumann and Dirichlet boundary macros
so that they can be used either for standard domain boundaries or for
embedded boundaries. The distinction between the two cases is based on
whether the `dirichlet` parameter is passed to the boundary function
(using the `data` parameter). */

@undef neumann
@def neumann(expr)   (data ? embed_area_center (point, &x, &y, &z),
		      *((bool *)data) = false, (expr) :
		      Delta*(expr) + val(_s,0,0,0))
@
@undef neumann_homogeneous
@def neumann_homogeneous() (data ? *((bool *)data) = false, (0) :
			    val(_s,0,0,0))
@
@undef dirichlet
@def dirichlet(expr) (data ? embed_area_center (point, &x, &y, &z),
		      *((bool *)data) = true, (expr) :
		      2.*(expr) - val(_s,0,0,0))
@
@undef dirichlet_homogeneous
@def dirichlet_homogeneous() (data ? *((bool *)data) = true, (0) :
			      - val(_s,0,0,0))
@

/**
## Prolongation for the multigrid solver

We use a simplified prolongation operator for the [multigrid
solver](poisson.h#mg_solve) i.e. simple injection if bilinear
interpolation would use values which are fully contained within the
embedded boundary. */

#if MULTIGRID
static inline double bilinear_embed (Point point, scalar s)
{
  if (!coarse(cs) || !coarse(cs,child.x))
    return coarse(s);
  #if dimension >= 2
  if (!coarse(cs,0,child.y) || !coarse(cs,child.x,child.y))
    return coarse(s);
  #endif
  #if dimension >= 3
  if (!coarse(cs,0,0,child.z) || !coarse(cs,child.x,0,child.z) ||
      !coarse(cs,0,child.y,child.z) ||
      !coarse(cs,child.x,child.y,child.z))
    return coarse(s);  
  #endif
  return bilinear (point, s);
}

#define bilinear(point, s) bilinear_embed(point, s)
#endif // MULTIGRID

/**
## Lifting the "small cell" CFL restriction

For explicit advection schemes, the timestep is limited by the CFL
conditions 
$$ 
\Delta t < \frac{c_s\Delta}{f_i|u_i|}
$$ 
where $i$ is the index of each face, and $c_s$ and $f_i$ are the
embedded volume and face fractions respectively. It is clear that the
timestep may need to be arbitrarily small if $c_s/f_s$ tends toward
zero. This is the "small cell" restriction of cut-cell finite-volume
techniques.

A classical technique to avoid this limitation is to use a "cell
merging" procedure, where the fluxes from cells which would "overflow"
are redistributed to neighboring cells.

The function below uses this approach to update a field *f*, advected
by the face velocity field *uf*, with corresponding advection fluxes
*flux*, during timestep *dt* which only verifies the standard CFL
condition
$$
\Delta t < \frac{\Delta}{|u_i|}
$$
*/

trace
void update_tracer (scalar f, face vector uf, face vector flux, double dt)
{

  /**
  Note that the distinction should be made between $c_m$, the cell
  fraction metric, and $c_s$, the embedded fraction. This is not done
  now so that embedded boundaries cannot be combined with a metric
  yet. 

  The field *e* will store the "overflowing" sum of fluxes for each cell. */
  
  scalar e[];
  foreach() {

    /**
    If the cell is empty, it cannot overflow. */
    
    if (cs[] <= 0.)
      e[] = 0.;

    /**
    If the cell does not contain an embedded boundary, it cannot
    overflow either and the sum of the fluxes can be added to advance
    *f* in time. */
    
    else if (cs[] >= 1.) {
      foreach_dimension()
	f[] += dt*(flux.x[] - flux.x[1])/Delta;
      e[] = 0.;
    }

    /**
    If the cell contains the embedded boundary, we compute the maximum
    timestep verifying the restrictive CFL condition
    $$ 
    \Delta t_{max} = \frac{c_s\Delta}{max(f_i|u_i|)}
    $$
    Note that *fs* does not appear in the code below because *uf*
    already stores the product $f_su$. */
    
    else {
      double umax = 0.;
      for (int i = 0; i <= 1; i++)
	foreach_dimension()
	  if (fabs(uf.x[i]) > umax)
	    umax = fabs(uf.x[i]);
      double dtmax = Delta*cs[]/(umax + SEPS);

      /**
      We compute the sum of the fluxes. */
      
      double F = 0.;
      foreach_dimension()
	F += flux.x[] - flux.x[1];
      F /= Delta*cs[];

      /**
      If the timestep is smaller than $\Delta t_{max}$, the cell
      cannot overflow and *f* is advanced in time using the entire
      flux. */
      
      if (dt <= dtmax) {
	f[] += dt*F;
	e[] = 0.;
      }
      
      /**
      Otherwise, the cell is filled "to the brim" by advancing *f*
      using the maximum allowable timestep. The field *e* is used to
      store the excess flux, weighted by the sum of the neighboring
      embedded fractions. */
      
      else {
	f[] += dtmax*F;
	double scs = 0.;
	foreach_neighbor(1)
	  scs += sq(cs[]);
	e[] = (dt - dtmax)*F*cs[]/scs;
      }
    }
  }
  boundary ({e});

  /**
  In a second phase, the excesses in each cell are added to the
  neighboring cells in proportion of their embedded fractions. */ 

  foreach() {
    double se = 0.;
    foreach_neighbor(1)
      se += e[];
    f[] += cs[]*se;
  }
}

/**
## Default settings

To apply the volume/area fraction-weighting to the solvers, we define
the metric using the embedded fractions. */

event metric (i = 0)
{
  foreach()
    cs[] = 1.;
  foreach_face()
    fs.x[] = 1.;
#if TREE
  cs.refine = embed_fraction_refine;

  /**
  For prolongation we cannot use the same function since the surface
  fraction field *fs* is not necessarily defined for prolongation
  cells. So we switch back to the default fraction refinement (which
  is less accurate but only relies on *cs*). */

  cs.prolongation = fraction_refine;
  foreach_dimension()
    fs.x.prolongation = embed_face_fraction_refine_x;
  
  /**
  Note that we do not need to change the `refine` method since the
  default `refine` method calls the prolongation method for each
  component. */
  
#endif
  boundary ({cs, fs});
  restriction ({cs, fs});

  // fixme: embedded boundaries cannot be combined with (another) metric yet
  assert (is_constant (cm) || cm.i == cs.i);
  
  cm = cs;
  fm = fs;
}

/**
## References

~~~bib
@article{johansen1998,
  title={A Cartesian grid embedded boundary method for Poisson's
  equation on irregular domains},
  author={Johansen, Hans and Colella, Phillip},
  journal={Journal of Computational Physics},
  volume={147},
  number={1},
  pages={60--85},
  year={1998},
  publisher={Elsevier},
  url={https://pdfs.semanticscholar.org/17cd/babecd054d58da05c2ba009cccb3c687f58f.pdf}
}

@article{schwartz2006,
  title={A Cartesian grid embedded boundary method for the heat equation 
  and Poisson’s equation in three dimensions},
  author={Schwartz, Peter and Barad, Michael and Colella, Phillip and Ligocki, 
  Terry},
  journal={Journal of Computational Physics},
  volume={211},
  number={2},
  pages={531--550},
  year={2006},
  publisher={Elsevier},
  url={https://cloudfront.escholarship.org/dist/prd/content/qt0fp606kk/qt0fp606kk.pdf}
}
~~~

## See also

* [Notes on drag force computation](/src/notes/drag.tm)
*/
