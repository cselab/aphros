#include "Cubism.h"

#include "../../hydro/vect.hpp"
#include "../../hydro/mesh3d.hpp"
#include "Hydro.h"

template <class KF>
Distr* Try(
    MPI_Comm comm, KernelFactory& kf, 
    int bs, Idx b, Idx p, int es, int h) {

  if (KF* r = dynamic_cast<KF*>(&kf)) {
    return new Cubism<KF>(comm, kf, bs, b, p, es, h);
  }
  return nullptr;
}

std::unique_ptr<Distr> CreateCubism(
    MPI_Comm comm, KernelFactory& kf, 
    int bs, Idx b, Idx p, int es, int h) {
  Distr* r = nullptr;
  //if (!r) r = Try<HydroFactory<geom::geom3d::MeshStructured<double>>(
  //    comm, kf, bs, b, p, es, h);
  //if (!r) r = Try<HydroFactory<geom::geom3d::MeshStructured<float>>(
  //    comm, kf, bs, b, p, es, h);
  assert(r && "CreateCubismz(): KernelFactory not found");
  return unique_ptr<Distr>(r);
}
