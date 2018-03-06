#include "Local.h"

#include "hydro/vect.hpp"
#include "hydro/mesh3d.hpp"
#include "Kernel.h"
#include "KernelMesh.h"
#include "Vars.h"

template <class KF>
Distr* TryLocal(
    MPI_Comm comm, KernelFactory& kf, int bs, int es, int h, Vars& par) {

  if (KF* kfd = dynamic_cast<KF*>(&kf)) {
    return new Local<KF>(comm, *kfd, bs, es, h, par);
  }
  return nullptr;
}

Distr* CreateLocal(
    MPI_Comm comm, KernelFactory& kf, int bs, int es, int h, Vars& par) {
  Distr* r = nullptr;
  if (!r) r = TryLocal<KernelMeshFactory<geom::MeshStructured<double, 3>>>(
      comm, kf, bs, es, h, par);
  //if (!r) r = Try<KernelMeshFactory<geom::geom3d::MeshStructured<float>>(
  //    comm, kf, bs, b, p, es, h);
  assert(r && "CreateLocal(): KernelFactory not found");
  return r;
}
