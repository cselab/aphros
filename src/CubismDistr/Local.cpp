#include "Local.h"

#include "hydro/vect.hpp"
#include "hydro/mesh3d.hpp"
#include "Kernel.h"
#include "KernelMesh.h"
#include "Vars.h"

template <class KF>
Distr* TryLocal(
    MPI_Comm comm, KernelFactory& kf, Vars& par) {

  if (KF* kfd = dynamic_cast<KF*>(&kf)) {
    return new Local<KF>(comm, *kfd, par);
  }
  return nullptr;
}

Distr* CreateLocal(
    MPI_Comm comm, KernelFactory& kf, Vars& par) {
  Distr* r = nullptr;
  if (!r) r = TryLocal<KernelMeshFactory<geom::MeshStructured<double, 3>>>(
      comm, kf, par);
  //if (!r) r = Try<KernelMeshFactory<geom::geom3d::MeshStructured<float>>(
  //    comm, kf, par);
  assert(r && "CreateLocal(): KernelFactory not found");
  return r;
}
