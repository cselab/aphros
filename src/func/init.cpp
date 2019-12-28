#include "init.h"
#include "geom/mesh.h"
#include "init_cl.h"
#include "init_sig.h"
#include "init_u.h"

using M = MeshStructured<double, 3>;

template void InitVf(FieldCell<typename M::Scal>& fcu, const Vars& var, M& m);

template std::function<void(
    FieldCell<typename M::Scal>&, const FieldCell<typename M::Scal>&, const M&)>
CreateInitCl(const Vars& par, bool verb);

template std::function<void(FieldCell<typename M::Scal>&, const M&)>
CreateInitU(const Vars& par, bool verb);

template std::function<void(FieldCell<typename M::Scal>&, const M&)>
CreateInitSig(const Vars& var);
