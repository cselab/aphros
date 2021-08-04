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
