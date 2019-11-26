// vim: expandtab:smarttab:sw=2:ts=2
#include <mpi.h>
#include <sstream>

#include "cubismnc.ipp"
#include "geom/mesh.h"
#include "kernel/kernel.h"
#include "kernel/kernelmesh.h"
#include "parse/vars.h"

template <size_t bx, size_t by, size_t bz, class KF>
Distr* Try(MPI_Comm comm, KernelFactory& kf, Vars& var) {
  using M = typename KF::M;
  using Scal = typename M::Scal;
  //using Par0 = GPar<Scal, bx, by, bz, 0>; // 0 halos on either side
  //using Par1 = GPar<Scal, bx, by, bz, 1>; // 1 halos on either side
  using Par2 = GPar<Scal, bx, by, bz, 2>; // 2 halos on either side
  //using Par3 = GPar<Scal, bx, by, bz, 3>; // 3 halos on either side

  // Check block size
  if (var.Int["bsx"] == bx &&
      var.Int["bsy"] == by &&
      var.Int["bsz"] == bz) {
    // Check kernel
    if (KF* kfd = dynamic_cast<KF*>(&kf)) {
      switch (var.Int["hl"]) {
        //case 0:
        //  return new Cubismnc<Par0, KF>(comm, *kfd, var);
        //case 1:
        //  return new Cubismnc<Par1, KF>(comm, *kfd, var);
        case 2:
          return new Cubismnc<Par2, KF>(comm, *kfd, var);
        //case 3:
        //  return new Cubismnc<Par3, KF>(comm, *kfd, var);
        default:
          break;
      }
    }
  }

  return nullptr;
}

Distr* CreateCubismnc(MPI_Comm comm, KernelFactory& kf, Vars& var) {
  Distr* r = nullptr;
  using KF = KernelMeshFactory<MeshStructured<double, 3>>;
  // 3D
  if (!r) r = Try<32, 32, 32, KF>(comm, kf, var);
  if (!r) r = Try<16, 16, 16, KF>(comm, kf, var);
  if (!r) r = Try<8, 8, 8, KF>(comm, kf, var);
  // 2D
  if (!r) r = Try<32, 32, 1, KF>(comm, kf, var);
  if (!r) r = Try<16, 16, 1, KF>(comm, kf, var);
  if (!r) r = Try<8, 8, 1, KF>(comm, kf, var);
  if (!r) {
    std::stringstream ss;
    ss << __func__ << ": no instance with "
      << "bs="
      << var.Int["bsx"] << " " << var.Int["bsy"] << " " << var.Int["bsz"]
      << " hl=" << var.Int["hl"];
    throw std::runtime_error(ss.str());
  }
  return r;
}

