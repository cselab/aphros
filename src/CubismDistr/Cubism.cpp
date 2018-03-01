#include "Cubism.h"

#include "hydro/vect.hpp"
#include "hydro/mesh3d.hpp"
#include "Kernel.h"
#include "KernelMesh.h"
#include "Vars.h"

template <class KF>
Distr* TryCubism(
    MPI_Comm comm, KernelFactory& kf, 
    int bs, int es, int h, Vars& par) {

  if (KF* kfd = dynamic_cast<KF*>(&kf)) {
    return new Cubism<KF>(comm, *kfd, bs, es, h, par);
  }
  return nullptr;
}

std::unique_ptr<Distr> CreateCubism(
    MPI_Comm comm, KernelFactory& kf, int bs, int es, int h, Vars& par) {
  Distr* r = nullptr;
  if (!r) r = TryCubism<KernelMeshFactory<geom::MeshStructured<double, 3>>>(
      comm, kf, bs, es, h, par);
  //if (!r) r = Try<KernelMeshFactory<geom::geom3d::MeshStructured<float>>(
  //    comm, kf, bs, b, p, es, h);
  assert(r && "CreateCubism(): KernelFactory not found");
  return unique_ptr<Distr>(r);
}
