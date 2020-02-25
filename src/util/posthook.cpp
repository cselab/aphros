#include "posthook.h"
#include "geom/mesh.h"

template <class M>
void PostHook(M&) {}

using M = MeshStructured<double, 3>;

template void PostHook(M& m);
