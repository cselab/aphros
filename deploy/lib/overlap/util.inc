static double
util_max(double x, double y)
{
  return x > y ? x : y;
}

static double
util_stable_norm(double x, double y, double z)
{
  double i, j;
  x = fabs(x);
  y = fabs(y);
  z = fabs(z);
  i = util_max(x, util_max(y, z));
  if (i == 0.0)
    return 0.0;
  else {
    j = 1/i;
    x *= j;
    y *= j;
    z *= j;
    return i * sqrt(x*x + y*y + z*z);
  }
}

static double
util_norm(double x, double y, double z)
{
  return sqrt(x*x + y*y + z*z);
}

static double
util_squared_norm(double x, double y, double z)
{
  return x*x + y*y + z*z;
}

/* s . ((a - b) ~ (x - y)) */
static double
util_diff_dot_cross(vector_t s,
		    vector_t a, vector_t b,
		    vector_t x, vector_t y) {
  double u, v, w;
  double h, j, k;

  u = a[0] - b[0];
  v = a[1] - b[1];
  w = a[2] - b[2];
  h = x[0] - y[0];
  j = x[1] - y[1];
  k = x[2] - y[2];
  return -s[0]*(j*w-k*v)-s[1]*(k*u-h*w)-s[2]*(h*v-j*u);
}
