#include <iostream>

#include <util/posthook.h>
#include <geom/mesh.h>

template <class M>
void PostHook(M&) {
  std::cout << "posthook_user_defined" << std::endl;
}

using M = MeshStructured<double, 3>;

template void PostHook(M& m);
