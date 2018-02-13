#include "Local.h"

#include "../../hydro/vect.hpp"
#include "../../hydro/mesh3d.hpp"
#include "Hydro.h"

template <class KF>
Distr* Try(KernelFactory& kf, int bs, Idx b, Idx p, int es, int h) {

  if (KF* kfd = dynamic_cast<KF*>(&kf)) {
    return new Local<KF>(*kfd, bs, b, p, es, h);
  }
  return nullptr;
}

std::unique_ptr<Distr> CreateLocal(
    KernelFactory& kf, int bs, Idx b, Idx p, int es, int h) {
  Distr* r = nullptr;
  r || (r = Try<HydroFactory<geom::geom3d::MeshStructured<double>>>(
      kf, bs, b, p, es, h));
  //if (!r) r = Try<HydroFactory<geom::geom3d::MeshStructured<float>>(
  //    comm, kf, bs, b, p, es, h);
  assert(r && "CreateLocalz(): KernelFactory not found");
  return std::unique_ptr<Distr>(r);
}
