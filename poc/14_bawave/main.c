}

event out (t += DUMPDT ; t <= TMAX + DUMPDT) {
   FILE * fp;    
  fprintf(stderr, "dump i=%05d t=%g dt=%g \n", i, t ,dt);

  static int frame = 0;
  scalar * a = {u, p, f};

  char name[1000];

  sprintf(name, "u_%04d.vtk", frame);
  fp = fopen(name, "w");
  io(a, fp);

  sprintf(name, "u_%04d.ppm", frame);
  fp = fopen(name, "w");
  output_ppm (f, fp, min = 0, max = 1, n = 512);
  fclose(fp);  

  ++frame;
}

event statout (i += 10 ; t <= TMAX) {
  fprintf(stderr, "step i=%05d t=%g dt=%g \n", i, t ,dt);
}

static double
wave0(double x, double y)
{
  double a_ = ak/k_;
  double eta1 = a_*cos(k_*x);
  double alpa = 1./tanh(k_*h_);
  double eta2 = 1./4.*alpa*(3.*sq(alpa) - 1.)*sq(a_)*k_*cos(2.*k_*x);
  double eta3 = -3./8.*(cube(alpa)*alpa - 
			3.*sq(alpa) + 3.)*cube(a_)*sq(k_)*cos(k_*x) + 
    3./64.*(8.*cube(alpa)*cube(alpa) + 
	    (sq(alpa) - 1.)*(sq(alpa) - 1.))*cube(a_)*sq(k_)*cos(3.*k_*x);
  return eta1 + ak*eta2 + sq(ak)*eta3 - y;
}

static double
wave(double x, double y)
{
    return wave0(x - 0.5, y - 0.5);
}

static double
phi0(double x, double y)
{
    double alpa = 1./tanh(k_*h_);
    double a_ = ak/k_;
    double sgma = sqrt(g_*k_*tanh(k_*h_)*
		       (1. + k_*k_*a_*a_*(9./8.*(sq(alpa) - 1.)*
					  (sq(alpa) - 1.) + sq(alpa))));
    double A_ = a_*g_/sgma;
    double phi1 = A_*cosh(k_*(y + h_))/cosh(k_*h_)*sin(k_*x);
    double phi2 = 3.*ak*A_/(8.*alpa)*(sq(alpa) - 1.)*(sq(alpa) - 1.)*
	cosh(2.0*k_*(y + h_))*sin(2.0*k_*x)/cosh(2.0*k_*h_);
    double phi3 = 1./64.*(sq(alpa) - 1.)*(sq(alpa) + 3.)*
	(9.*sq(alpa) - 13.)*
	cosh(3.*k_*(y + h_))/cosh(3.*k_*h_)*a_*a_*k_*k_*A_*sin(3.*k_*x);
    return phi1 + ak*phi2 + ak*ak*phi3;
}

static double
phi(double x, double y)
{
    return phi0(x - 0.5, y - 0.5);
}

