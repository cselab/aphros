#include "init_u.h"
#include "geom/mesh.h"

using M = MeshStructured<double, 3>;

template void InitVf(FieldCell<typename M::Scal>& fcu, const Vars& var, M& m);
