#include "geom/mesh.h"
#include "solver/reconst.h"

using Scal = double;
using R = Reconst<Scal>;
using std::cout;

int main() {
    enum {X, Y, Z};
    Scal u, u0, a;
    GVect<Scal, 3> n;

    n[X] = 0.2; n[Y] = 0.3; n[Z] = 0;
    u = 0.1;
    a = R::GetLineA1(n, u);
    u0 = R::GetLineU1(n, a);
    cout << a << ' ' << u << ' ' << u0 << '\n';
    return 0;
}
