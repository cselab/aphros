#include "geom/vect.h"

// Returns volume of intersection of cell and sphere.
// x: cell center
// h: cell size
// xc: sphere center
// r: radius
double GetSphereOverlap(const GVect<double, 3>& x, const GVect<double, 3>& h, 
                        const GVect<double, 3>& xc, double r);
