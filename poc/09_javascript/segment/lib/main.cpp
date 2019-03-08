#include <solver/normal.h>
#include <solver/reconst.h>
#include <geom/mesh.h>

#include "segment.h"

const int D = 5, dim = 3, edim = 2;
using MIdx = GMIdx<dim>;
using Dir = GDir<dim>;
using Vect = GVect<Scal, dim>;
using M = MeshStructured<Scal, dim>;
using R = Reconst<Scal>;
using U = solver::UNormal<M>;

struct Segment {
    Scal n[2*D*D], a[D*D];
} segment;

int segment_get(const Scal alpha[D*D], /**/ Scal **pn, Scal **pa) {
    enum {X, Y, Z};
    Rect<Vect> dom(Vect(0), Vect(1));
    MIdx b(0);
    MIdx s(D, D, 1);
    int hl;
    FieldCell<Vect> fcn;
    FieldCell<Scal> fck;
    Vect u;
    MIdx w;
    Scal *n, *a, al;
    int i;

    hl = 2;
    M m = InitUniformMesh<M>(dom, b, s, hl, true, s);
    FieldCell<Scal> fcu(m, 0);
    FieldCell<bool> fci(m, true);

    n = segment.n;
    a = segment.a;

    i = 0;
    for (auto c : m.Cells()) {
      fcu[c] = alpha[i++];
    }

    U::CalcNormal(m, fcu, fci, edim, /**/ fcn, fck);

    i = 0;
    for (auto c : m.Cells()) {
      u = fcn[c];
      n[i]       = u[X];
      n[i + D*D] = u[Y];
      i++;
    }

    i = 0;
    for (auto c : m.Cells()) {
        u = fcn[c];
        al = alpha[i];
        a[i] = R::GetLineA1(u, al);
        i++;
    }

    *pn = n; *pa = a;
    return 0;
}
