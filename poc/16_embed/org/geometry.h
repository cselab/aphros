/**
# Basic geometric functions

These basic geometric functions are mostly related to Volume-Of-Fluid
computations.

We consider a square cell of size unity centered on the origin, cut by
a straight line.

![Cell and interface](/src/figures/square.svg)

The line can be described by the equation
$$
n_xx+n_yy=\alpha
$$
where $\mathbf{n}$ is a vector normal to the interface and $\alpha$ is
the intercept. We note $c$ the volume of the part of the square cell
which lies "inside" the interface, where "inside" is defined by
convention as the opposite direction to the normal vector $\mathbf{n}$
(i.e. the normal vector is pointing "outside").

With these definitions, the interface is uniquely defined by providing
$\mathbf{n}$ and either $\alpha$ or $c$ i.e. there is a unique
function which computes $\alpha$ given $c$ and $\mathbf{n}$. We call
this function `line_alpha()` and define it as: */

#if dimension >= 2
double line_alpha (double c, coord n)
{
  double alpha, n1, n2;
  
  n1 = fabs (n.x); n2 = fabs (n.y);
  if (n1 > n2)
    swap (double, n1, n2);

  c = clamp (c, 0., 1.);
  double v1 = n1/2.;
  if (c <= v1/n2)
    alpha = sqrt (2.*c*n1*n2);
  else if (c <= 1. - v1/n2)
    alpha = c*n2 + v1;
  else
    alpha = n1 + n2 - sqrt (2.*n1*n2*(1. - c));

  if (n.x < 0.)
    alpha += n.x;
  if (n.y < 0.)
    alpha += n.y;

  return alpha - (n.x + n.y)/2.;
}
#endif // dimension >= 2

#if dimension >= 3
double plane_alpha (double c, coord n)
{
  double alpha;
  coord n1;
  
  n1.x = fabs (n.x); n1.y = fabs (n.y); n1.z = fabs (n.z);

  double m1, m2, m3;
  m1 = min(n1.x, n1.y);
  m3 = max(n1.x, n1.y);
  m2 = n1.z;
  if (m2 < m1) {
    double tmp = m1;
    m1 = m2;
    m2 = tmp;
  }
  else if (m2 > m3) {
    double tmp = m3;
    m3 = m2;
    m2 = tmp;
  }
  double m12 = m1 + m2;
  double pr = max(6.*m1*m2*m3, 1e-50);
  double V1 = m1*m1*m1/pr;
  double V2 = V1 + (m2 - m1)/(2.*m3), V3;
  double mm;
  if (m3 < m12) {
    mm = m3;
    V3 = (m3*m3*(3.*m12 - m3) + m1*m1*(m1 - 3.*m3) + m2*m2*(m2 - 3.*m3))/pr;
  }
  else {
    mm = m12;
    V3 = mm/(2.*m3);
  }

  c = clamp (c, 0., 1.);
  double ch = min(c, 1. - c);
  if (ch < V1)
    alpha = pow (pr*ch, 1./3.);
  else if (ch < V2)
    alpha = (m1 + sqrt(m1*m1 + 8.*m2*m3*(ch - V1)))/2.;
  else if (ch < V3) {
    double p12 = sqrt (2.*m1*m2);
    double q = 3.*(m12 - 2.*m3*ch)/(4.*p12);
    double teta = acos(clamp(q,-1.,1.))/3.;
    double cs = cos(teta);
    alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + m12;
  }
  else if (m12 < m3)
    alpha = m3*ch + mm/2.;
  else {
    double p = m1*(m2 + m3) + m2*m3 - 1./4., p12 = sqrt(p);
    double q = 3.*m1*m2*m3*(1./2. - ch)/(2.*p*p12);
    double teta = acos(clamp(q,-1.,1.))/3.;
    double cs = cos(teta);
    alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + 1./2.;
  }
  if (c > 1./2.) alpha = 1. - alpha;

  if (n.x < 0.)
    alpha += n.x;
  if (n.y < 0.)
    alpha += n.y;
  if (n.z < 0.)
    alpha += n.z;

  return alpha - (n.x + n.y + n.z)/2.;;
}
#else // dimension < 3
# define plane_alpha line_alpha
#endif

/**
Conversely there is a unique function computing $c$ as a function of
$\mathbf{n}$ and $\alpha$. We call this function `line_area()` and
define it as: */

#if dimension >= 2
double line_area (double nx, double ny, double alpha)
{
  double a, v, area;

  alpha += (nx + ny)/2.;
  if (nx < 0.) {
    alpha -= nx;
    nx = - nx;
  }
  if (ny < 0.) {
    alpha -= ny;
    ny = - ny;
  }

  if (alpha <= 0.)
    return 0.;

  if (alpha >= nx + ny)
    return 1.;

  if (nx < 1e-10)
    area = alpha/ny;
  else if (ny < 1e-10)
    area = alpha/nx;
  else {
    v = sq(alpha);

    a = alpha - nx;
    if (a > 0.)
      v -= a*a;
    
    a = alpha - ny;
    if (a > 0.)
      v -= a*a;

    area = v/(2.*nx*ny);
  }

  return clamp (area, 0., 1.);
}
#endif // dimension >= 2

#if dimension >= 3
double plane_volume (coord n, double alpha)
{
  double al = alpha + (n.x + n.y + n.z)/2. +
    max(0., -n.x) + max(0., -n.y) + max(0., -n.z);
  if (al <= 0.)
    return 0.;
  double tmp = fabs(n.x) + fabs(n.y) + fabs(n.z);
  if (al >= tmp)
    return 1.;
  if (tmp < 1e-10)
    return 0.;
  double n1 = fabs(n.x)/tmp;
  double n2 = fabs(n.y)/tmp;
  double n3 = fabs(n.z)/tmp;
  al = max(0., min(1., al/tmp));
  double al0 = min(al, 1. - al);
  double b1 = min(n1, n2);
  double b3 = max(n1, n2);
  double b2 = n3;
  if (b2 < b1) {
    tmp = b1;
    b1 = b2;
    b2 = tmp;
  }
  else if (b2 > b3) {
    tmp = b3;
    b3 = b2;
    b2 = tmp;
  }
  double b12 = b1 + b2;
  double bm = min(b12, b3);
  double pr = max(6.*b1*b2*b3, 1e-50);
  if (al0 < b1)
    tmp = al0*al0*al0/pr;
  else if (al0 < b2)
    tmp = 0.5*al0*(al0 - b1)/(b2*b3) +  b1*b1*b1/pr;
  else if (al0 < bm)
    tmp = (al0*al0*(3.*b12 - al0) + b1*b1*(b1 - 3.*al0) +
	   b2*b2*(b2 - 3.*al0))/pr;
  else if (b12 < b3)
    tmp = (al0 - 0.5*bm)/b3;
  else
    tmp = (al0*al0*(3. - 2.*al0) + b1*b1*(b1 - 3.*al0) + 
	   b2*b2*(b2 - 3.*al0) + b3*b3*(b3 - 3.*al0))/pr;

  double volume = al <= 0.5 ? tmp : 1. - tmp;
  return clamp (volume, 0., 1.);
}
#else // dimension < 3
# define plane_volume(n, alpha) line_area(n.x, n.y, alpha)
#endif

/**
VOF algorithms require the computation of volume fractions on
(rectangular) parts of the initial square cell.

We first define a function which takes an interface definition
($\mathbf{n}$, $\alpha$), the coordinates of the lower-left `a`
and upper-right `b` corners of a rectangle and returns the
fraction of this rectangle which lies inside the interface. */

double rectangle_fraction (coord n, double alpha, coord a, coord b)
{
  coord n1;
  foreach_dimension() {
    alpha -= n.x*(b.x + a.x)/2.;
    n1.x = n.x*(b.x - a.x);
  }
  return plane_volume (n1, alpha);
}

/**
From the interface definition, it is also possible to compute the
coordinates of the segment in 2D, or facet in 3D, representing the
interface in the unit cell.

In two dimensions, the function below returns the 0,1 or 2 coordinates
(stored in the `p` array provided by the user) of the corresponding
interface segments. The case where only 1 coordinate is returned
corresponds to the degenerate case where the interface intersects the
cell exactly on a vertex. 

In three dimensions, the function returns up to 12 coordinates of the
planar fragment. */

#if dimension <= 2
int facets (coord n, double alpha, coord p[2])
{
  int i = 0;
  for (double s = -0.5; s <= 0.5; s += 1.)
    foreach_dimension()
      if (fabs (n.y) > 1e-4 && i < 2) {
	double a = (alpha - s*n.x)/n.y;
	if (a >= -0.5 && a <= 0.5) {
	  p[i].x   = s;
	  p[i++].y = a;
	}
      }
  return i;
}
#else
static coord cube_edge[12][2] = {
  {{0.,0.,0.},{1.,0.,0.}},{{0.,0.,1.},{1.,0.,1.}},
  {{0.,1.,1.},{1.,1.,1.}},{{0.,1.,0.},{1.,1.,0.}},
  {{0.,0.,0.},{0.,1.,0.}},{{0.,0.,1.},{0.,1.,1.}},
  {{1.,0.,1.},{1.,1.,1.}},{{1.,0.,0.},{1.,1.,0.}},
  {{0.,0.,0.},{0.,0.,1.}},{{1.,0.,0.},{1.,0.,1.}},
  {{1.,1.,0.},{1.,1.,1.}},{{0.,1.,0.},{0.,1.,1.}}
};

/* first index is the edge number, second index is the edge
   orientation (0 or 1), third index are the edges which this edge may
   connect to in order */
static int cube_connect[12][2][4] = {
  {{9, 1, 8}, {4, 3, 7}},   /* 0 */
  {{6, 2, 5},  {8, 0, 9}}, /* 1 */
  {{10, 3, 11},  {5, 1, 6}},  /* 2 */
  {{7, 0, 4},   {11, 2, 10}},  /* 3 */
  {{3, 7, 0},   {8, 5, 11}},  /* 4 */
  {{11, 4, 8},  {1, 6, 2}},  /* 5 */
  {{2, 5, 1},  {9, 7, 10}}, /* 6 */
  {{10, 6, 9}, {0, 4, 3}},   /* 7 */
  {{5, 11, 4},  {0, 9, 1}}, /* 8 */
  {{1, 8, 0}, {7, 10, 6}}, /* 9 */
  {{6, 9, 7},  {3, 11, 2}},   /* 10 */
  {{2, 10, 3},   {4, 8, 5}}    /* 11 */
};

int facets (coord n, double alpha, coord v[12], double h)
{
  coord a[12];
  int orient[12];

  for (int i = 0; i < 12; i++) {
    coord e, d;
    double den = 0., t = alpha;
    foreach_dimension() {
      d.x = h*(cube_edge[i][0].x - 0.5);
      e.x = h*(cube_edge[i][1].x - 0.5);
      den += n.x*(e.x - d.x);
      t -= n.x*d.x;
    }
    orient[i] = -1;
    if (fabs (den) > 1e-10) {
      t /= den;
      if (t >= 0. && t < 1.) {
	double s = - alpha;
	foreach_dimension() {
	  a[i].x = d.x + t*(e.x - d.x);
	  s += n.x*e.x;
	}
	orient[i] = (s > 0.);
      }
    }
  }

  for (int i = 0; i < 12; i++) {
    int nv = 0, e = i;
    while (orient[e] >= 0) {
      int m = 0, * ne = cube_connect[e][orient[e]];
      v[nv++] = a[e];
      orient[e] = -1;
      while (m < 3 && orient[e] < 0)
	e = ne[m++];
    }
    if (nv > 2)
      return nv;
  }
  return 0;
}
#endif // dimension == 3

/**
This function fills the coordinates *p* of the centroid of the
interface fragment and returns the length/area of the fragment. */

double line_length_center (coord m, double alpha, coord * p)
{
  alpha += (m.x + m.y)/2.;

  coord n = m;
  foreach_dimension(2)
    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }

  p->x = p->y = p->z = 0.;

  if (alpha <= 0. || alpha >= n.x + n.y)
    return 0.;
  
  foreach_dimension(2)
    if (n.x < 1e-4) {
      p->x = 0.;
      p->y = (m.y < 0. ? 1. - alpha : alpha) - 0.5;
      return 1.;
    }
  
  if (alpha >= n.x) {
    p->x += 1.;
    p->y += (alpha - n.x)/n.y;
  }
  else
    p->x += alpha/n.x;

  double ax = p->x, ay = p->y;
  if (alpha >= n.y) {
    p->y += 1.;
    ay -= 1.;
    p->x += (alpha - n.y)/n.x;
    ax -= (alpha - n.y)/n.x;
  }
  else {
    p->y += alpha/n.y;
    ay -= alpha/n.y;
  }

  foreach_dimension(2) {
    p->x /= 2.;
    p->x = clamp (p->x, 0., 1.);
    if (m.x < 0.)
      p->x = 1. - p->x;
    p->x -= 0.5;
  }

  return sqrt (ax*ax + ay*ay);
}

#if dimension == 2
#define plane_area_center(m,a,p) line_length_center(m,a,p)
#else // dimension == 3
double plane_area_center (coord m, double alpha, coord * p)
{
  foreach_dimension()
    if (fabs (m.x) < 1e-4) {
      coord n, q;
      ((double *)&n)[0] = m.y;
      ((double *)&n)[1] = m.z;
      double length = line_length_center (n, alpha, &q);
      p->x = 0.;
      p->y = ((double *)&q)[0];
      p->z = ((double *)&q)[1];
      return length;
    }

  alpha += (m.x + m.y + m.z)/2.;
  coord n = m;
  foreach_dimension()
    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }

  double amax = n.x + n.y + n.z;
  if (alpha < 0. || alpha > amax) {
    p->x = p->y = p->z = 0.;
    return 0.;
  }

  double area = sq(alpha);
  p->x = p->y = p->z = area*alpha;

  foreach_dimension() {
    double b = alpha - n.x;
    if (b > 0.) {
      area -= b*b;
      p->x -= b*b*(2.*n.x + alpha);
      p->y -= b*b*b;
      p->z -= b*b*b;
    }
  }

  amax = alpha - amax;
  foreach_dimension() {
    double b = amax + n.x;
    if (b > 0.) {
      area += b*b;
      p->y += b*b*(2.*n.y + alpha - n.z);
      p->z += b*b*(2.*n.z + alpha - n.y);
      p->x += b*b*b;
    }
  }

  area *= 3.;
  foreach_dimension() {
    if (area) {
      p->x /= area*n.x;
      p->x = clamp (p->x, 0., 1.);
    }
    else
      p->x = 0.;
    if (m.x < 0.) p->x = 1. - p->x;
    p->x -= 0.5;
  }

  return area*sqrt (1./(sq(n.x)*sq(n.y)) +
		    1./(sq(n.x)*sq(n.z)) +
		    1./(sq(n.z)*sq(n.y)))/6.;
}
#endif // dimension == 3

/**
This function fills the coordinates *p* of the centroid of the
fraction *a* of a square cell lying under the line $(m,alpha)$. */

void line_center (coord m, double alpha, double a, coord * p)
{
  alpha += (m.x + m.y)/2.;
  
  coord n = m;
  foreach_dimension(2)
    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }

  p->z = 0.;
  if (alpha <= 0.) {
    p->x = p->y = -0.5;
    return;
  }

  if (alpha >= n.x + n.y) {
    p->x = p->y = 0.;
    return;
  }

  foreach_dimension(2)
    if (n.x < 1e-4) {
      p->x = 0.;
      p->y = sign(m.y)*(a/2. - 0.5);
      return;
    }

  p->x = p->y = cube(alpha);
  
  foreach_dimension(2) {
    double b = alpha - n.x;
    if (b > 0.) {
      p->x -= sq(b)*(alpha + 2.*n.x);
      p->y -= cube(b);
    }
  }

  foreach_dimension(2) {
    p->x /= 6.*sq(n.x)*n.y*a;
    p->x = sign(m.x)*(p->x - 0.5);
  }
}

/**
This function fills the coordinates *p* of the centroid of the
fraction *a* of a cubic cell lying under the plane $(m,alpha)$. */

#if dimension == 2
#define plane_center(m,alpha,a,p) line_center(m,alpha,a,p)
#else // dimension == 3
void plane_center (coord m, double alpha, double a, coord * p)
{
  foreach_dimension()
    if (fabs (m.x) < 1e-4) {
      coord n, q;
      ((double *)&n)[0] = m.y;
      ((double *)&n)[1] = m.z;
      line_center (n, alpha, a, &q);
      p->x = 0.;
      p->y = ((double *)&q)[0];
      p->z = ((double *)&q)[1];
      return;
    }

  alpha += (m.x + m.y + m.z)/2.;
  coord n = m;
  foreach_dimension()
    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }

  if (alpha <= 0. || a == 0.) {
    p->x = p->y = p->z = -0.5;
    return;
  }

  if (alpha >= n.x + n.y + n.z || a == 1.) {
    p->x = p->y = p->z = 0.;
    return;
  }

  p->x = p->y = p->z = sq(sq(alpha));
  foreach_dimension() {
    double b = alpha - n.x;
    if (b > 0.) {
      p->x -= cube(b)*(3.*n.x + alpha);
      p->y -= sq(sq(b));
      p->z -= sq(sq(b));
    }
  }

  double amax = alpha - (n.x + n.y + n.z);
  foreach_dimension() {
    double b = amax + n.z;
    if (b > 0.) {
      p->x += cube(b)*(3.*n.x + alpha - n.y);
      p->y += cube(b)*(3.*n.y + alpha - n.x);
      p->z += sq(sq(b));
    }
  }

  double b = 24.*n.x*n.y*n.z*a;
  foreach_dimension() {
    p->x /= b*n.x;
    p->x = sign(m.x)*(p->x - 0.5);
  }
}
#endif // dimension == 3
