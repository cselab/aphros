#include "Cubism.h"

#include "hydro/vect.hpp"
#include "hydro/mesh3d.hpp"
#include "Kernel.h"
#include "KernelMesh.h"
#include "Vars.h"

template <class KF>
Distr* TryCubism(
    MPI_Comm comm, KernelFactory& kf, Vars& par) {

  if (KF* kfd = dynamic_cast<KF*>(&kf)) {
    return new Cubism<KF>(comm, *kfd, par);
  }
  return nullptr;
}

Distr* CreateCubism(
    MPI_Comm comm, KernelFactory& kf, Vars& par) {
  Distr* r = nullptr;
  if (!r) r = TryCubism<KernelMeshFactory<geom::MeshStructured<double, 3>>>(
      comm, kf, par);
  //if (!r) r = Try<KernelMeshFactory<geom::geom3d::MeshStructured<float>>(
  //    comm, kf, par);
  assert(r && "CreateCubism(): KernelFactory not found");
  return r;
}

std::string StreamHdfDyn::NAME = "alpha";
const std::string StreamHdfDyn::EXT = "";
int StreamHdfDyn::ID = 0;
