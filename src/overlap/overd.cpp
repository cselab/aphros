#include "over.h"

// Dummy implementation 

double GetSphereOverlap(const GVect<double, 3>& x, const GVect<double, 3>& h, 
                        const GVect<double, 3>& c, double r) {
  double d = (r - x.dist(c)) / h.norm();
  d = (d + 1.) * 0.5;
  d = std::min(d, 1.);
  d = std::max(d, 0.);
  return d;
}
