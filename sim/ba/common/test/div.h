#define SA 3
#define SW 1
#define DIM (SA*SA)
#define GET(x,y) (u[(y)*SA+(x)])
#define U(x,y) (u[(y+SW)*SA+(x+SW)])

static double Norm(double x, double y) {
  return sqrt(x * x + y * y);
}

static double Curv0(double u[]) {
  // gradient
  double gxpx = U(1,0) - U(0,0);
  double gxpy = ((U(0,1) + U(1,1)) - (U(0,-1) + U(1,-1))) * 0.25;

  double gxmx = U(0,0) - U(-1,0);
  double gxmy = ((U(-1,1) + U(0,1)) - (U(-1,-1) + U(0,-1))) * 0.25;

  double gypx = ((U(1,0) + U(1,1)) - (U(-1,0) + U(-1,1))) * 0.25;
  double gypy = U(0,1) - U(0,0);

  double gymx = ((U(1,-1) + U(1,0)) - (U(-1,-1) + U(-1,0))) * 0.25;
  double gymy = U(0,0) - U(0,-1);

  // length
  double gxp = Norm(gxpx, gxpy);
  double gxm = Norm(gxmx, gxmy);
  double gyp = Norm(gypx, gypy);
  double gym = Norm(gymx, gymy);

  // divergence of unit normal
  double d = gxpx / gxp - gxmx / gxm + gypy / gyp - gymy / gym;
  return d;
}
