#include "over.h"

double GetSphereOverlap(const GVect<double, 3>& x, const GVect<double, 3>& h, 
                        const GVect<double, 3>& xc, double r) {
  return x[0] + h[0] + xc[0] + r + 1.;
}
