Scal Clip(Scal a, Scal l, Scal u) {
  return std::max(l, std::min(u, a));
}

/**
 * gfs_line_area:
 * @m: normal to the line.
 * @alpha: line constant.
 *
 * Returns: the area of the fraction of a cell lying under the line
 * (@m,@alpha).
 */
Scal gfs_line_area (const Vect& m, Scal alpha)
{
  Vect n;
  Scal alpha1, a, v, area;

  g_return_val_if_fail (m != NULL, 0.);

  n = *m;
  alpha1 = alpha;
  if (n.x < 0.) {
    alpha1 -= n.x;
    n.x = - n.x;
  }
  if (n.y < 0.) {
    alpha1 -= n.y;
    n.y = - n.y;
  }

  if (alpha1 <= 0.)
    return 0.;

  if (alpha1 >= n.x + n.y)
    return 1.;

  if (n.x == 0.)
    area = alpha1/n.y;
  else if (n.y == 0.)
    area = alpha1/n.x;
  else {
    v = alpha1*alpha1;

    a = alpha1 - n.x;
    if (a > 0.)
      v -= a*a;

    a = alpha1 - n.y;
    if (a > 0.)
      v -= a*a;

    area = v/(2.*n.x*n.y);
  }

  return Clip(area, 0., 1.);
}

/**
 * gfs_line_alpha:
 * @m: a #Vect.
 * @c: a volume fraction.
 *
 * Returns: the value @alpha such that the area of a square cell
 * lying under the line defined by @m.@x = @alpha is equal to @c. 
 */
Scal gfs_line_alpha (const Vect& m, Scal c)
{
  Scal alpha, m1, m2, v1;

  g_return_val_if_fail (m != NULL, 0.);
  g_return_val_if_fail (c >= 0. && c <= 1., 0.);

  m1 = fabs (m->x); m2 = fabs (m->y);
  if (m1 > m2) {
    v1 = m1; m1 = m2; m2 = v1;
  }

  v1 = m1/2.;
  if (c <= v1/m2)
    alpha = sqrt (2.*c*m1*m2);
  else if (c <= 1. - v1/m2)
    alpha = c*m2 + v1;
  else
    alpha = m1 + m2 - sqrt (2.*m1*m2*(1. - c));

  if (m->x < 0.)
    alpha += m->x;
  if (m->y < 0.)
    alpha += m->y;

  return alpha;
}

#define EPS 1e-4

/**
 * gfs_line_center:
 * @m: normal to the line.
 * @alpha: line constant.
 * @a: area of cell fraction.
 * @p: a #Vect.
 *
 * Fills @p with the position of the center of mass of the fraction of
 * a square cell lying under the line (@m,@alpha).
 */
void gfs_line_center (const Vect& m, Scal alpha, Scal a, Vect * p)
{
  Vect n;
  Scal b;

  g_return_if_fail (m != NULL);
  g_return_if_fail (p != NULL);

  n = *m;
  if (n.x < 0.) {
    alpha -= n.x;
    n.x = - n.x;
  }
  if (n.y < 0.) {
    alpha -= n.y;
    n.y = - n.y;
  }

  p->z = 0.;
  if (alpha <= 0.) {
    p->x = p->y = 0.;
    return;
  }

  if (alpha >= n.x + n.y) {
    p->x = p->y = 0.5;
    return;
  }

  g_return_if_fail (a > 0. && a < 1.);

  if (n.x < EPS) {
    p->x = 0.5;
    p->y = m->y < 0. ? 1. - a/2. : a/2.;
    return;
  }

  if (n.y < EPS) {
    p->y = 0.5;
    p->x = m->x < 0. ? 1. - a/2. : a/2.;
    return;
  }

  p->x = p->y = alpha*alpha*alpha;

  b = alpha - n.x;
  if (b > 0.) {
    p->x -= b*b*(alpha + 2.*n.x);
    p->y -= b*b*b;
  }

  b = alpha - n.y;
  if (b > 0.) {
    p->y -= b*b*(alpha + 2.*n.y);
    p->x -= b*b*b;
  }

  p->x /= 6.*n.x*n.x*n.y*a;
  p->y /= 6.*n.x*n.y*n.y*a;

  if (m->x < 0.)
    p->x = 1. - p->x;
  if (m->y < 0.)
    p->y = 1. - p->y;
}

/**
 * gfs_line_area_center:
 * @m: normal to the line.
 * @alpha: line constant.
 * @p: a #Vect.
 *
 * Fills @p with the position of the center of area of the fraction of
 * a square cell lying under the line (@m,@alpha).
 *
 * Returns: the length of the facet.
 */
Scal gfs_line_area_center (const Vect& m, Scal alpha, Vect * p)
{
  Vect n;

  g_return_val_if_fail (m != NULL, 0.);
  g_return_val_if_fail (p != NULL, 0.);

  n = *m;
  if (n.x < 0.) {
    alpha -= n.x;
    n.x = - n.x;
  }
  if (n.y < 0.) {
    alpha -= n.y;
    n.y = - n.y;
  }

  p->z = 0.;
  if (alpha <= 0. || alpha >= n.x + n.y) {
    p->x = p->y = 0.;
    return 0.;
  }

  if (n.x < EPS) {
    p->x = 0.5;
    p->y = m->y < 0. ? 1. - alpha : alpha;
    return 1.;
  }

  if (n.y < EPS) {
    p->y = 0.5;
    p->x = m->x < 0. ? 1. - alpha : alpha;
    return 1.;
  }

  p->x = p->y = 0.;

  if (alpha >= n.x) {
    p->x += 1.;
    p->y += (alpha - n.x)/n.y;
  }
  else
    p->x += alpha/n.x;

  Scal ax = p->x, ay = p->y;
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

  p->x /= 2.;
  p->y /= 2.;

  THRESHOLD (p->x);
  THRESHOLD (p->y);

  if (m->x < 0.)
    p->x = 1. - p->x;
  if (m->y < 0.)
    p->y = 1. - p->y;

  return sqrt (ax*ax + ay*ay);
}
